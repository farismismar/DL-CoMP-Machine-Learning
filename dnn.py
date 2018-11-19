# -*- coding: utf-8 -*-
"""
Created on Thu April 28 20:46:45 2017

@author: farismismar
"""
# # Use these lines to force GPU
# import os
# os.environ["CUDA_DEVICE_ORDER"]="PCI_BUS_ID"
# # The GPU id to use, usually either "0" or "1"
# os.environ["CUDA_VISIBLE_DEVICES"]="0" 
# ###################################################

# # Do other imports now...
# import tensorflow as tf
# config = tf.ConfigProto()
# config.gpu_options.per_process_gpu_memory_fraction = 0.3

# from keras.backend.tensorflow_backend import set_session
# set_session(tf.Session(config=config))

# import keras

from keras.models import Sequential
from keras.layers import Dense

from sklearn.model_selection import GridSearchCV
from keras.wrappers.scikit_learn import KerasClassifier
#import keras.backend as K

import pandas as pd
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
#from matplotlib import rc
import numpy as np
from sklearn.preprocessing import MinMaxScaler

from sklearn.metrics import roc_curve, auc

from scipy.io import loadmat

# Set the random seed
seed = 3
import tensorflow as tf
tf.set_random_seed(seed)

# save X and Y dimensions   
mX, nX = [0,0]
mY, nY = [0,0]

y_pred = None # placeholder

episode_count = 0
model = None # the forward definition.
ss = None # The scaler forward def

def initialize_wrapper(random_state):
    global seed
    seed = random_state
    np.random.seed(random_state)
    
def predict_wrapper(filename): # filename = newX.mat
    data = loadmat(filename) # this is a dict.

#    keys = list(data.keys())[3:] # skip the first three columns
    values = list(data.values())[3:]
    
    a = []
    for sublist in values:
        for items in sublist:
            a.append(items)
            
    X_test = pd.DataFrame(a)
    del a
    X_test.columns = ['SINR', 'RSRP']
    
    global model
    global y_pred
    global ss

    return model.predict(ss.transform(X_test))

# returns the test error
def train_wrapper(filename): # filename='measurements.mat'
    global dataset
    global seed
    global model
    global ss
    global mX
    global nX
    global mY
    global nY
    global episode_count
     
    data = loadmat(filename) # this is a dict.
    keys = list(data.keys())[3:] # skip the first three columns
    values = list(data.values())[3:]
    
    dataset = pd.DataFrame()
    dataset = dataset.reindex(columns = keys) # create an empty dataframe

    for ii in np.arange(len(values)):
        v_ = np.array(values[ii])
        dataset[keys[ii]] = pd.Series(v_.flatten()) # cannot add the data to this empty df.
        

#    print(dataset)
    
    dataset['y'] = 1*(dataset['BLER'] <= 0.1) # H-ARQ target.
    dataset = dataset.drop(['BLER', 'CQI'], axis=1)

    # Perform a split 30-70
    train, test = train_test_split(dataset, test_size = 0.30, random_state = seed)

    # TO DO
    # Ordinal variables
    # No need for dummy coding here, just a quick replace.
    #train['Status'].replace(['Show-Up', 'No-Show'], [0, 1],inplace=True)
    #test['Status'].replace(['Show-Up', 'No-Show'], [0, 1],inplace=True)
    
    X_train = train.drop('y', axis = 1)
    X_test = test.drop('y', axis = 1)
    
    y_train = train['y'].values
    y_test = test['y'].values

    mX, nX = X_train.shape
    mY = y_train.shape
    nY = 1
    
    # Now undersample
    try:
        train_c0 = train[(y_train == 0)]
        train_c1= train[(y_train == 1)].sample(train_c0.shape[0], random_state=seed)
        train_undersample = train_c1.append(train_c0)
        X_train = train_undersample.drop('y', axis = 1)
        y_train = train_undersample['y']
    except:
        print('Undersampling failed')
        return [None, None, -1] # model is invalid

    ss = MinMaxScaler(feature_range=(0,1))
    
    # Scale the variables
    X_train_sc = ss.fit_transform(X_train)
    X_test_sc = ss.transform(X_test)
    
    model = KerasClassifier(build_fn=create_mlp, verbose=0, epochs=5, batch_size=32)

    # The hyperparameters
    width_dims=[3,5]
    n_hiddens = [1,3] # the depth of hidden layers
    activators = ['sigmoid', 'relu'] # 'softmax']

    hyperparameters = dict(width=width_dims, depth=n_hiddens, act=activators)
    
    grid = GridSearchCV(estimator=model, param_grid=hyperparameters, n_jobs=1, cv=3)
    grid_result = grid.fit(X_train_sc, y_train)
    
    # This is the best model
    best_model_mlp = grid_result.best_params_
    print(best_model_mlp)
 
    mlp = grid_result.best_estimator_
    mlp.fit(X_train_sc, y_train)
    model = mlp # the final model
    
#    y_pred_dnn = mlp.predict(X_test_sc)
    y_score_dnn = mlp.predict_proba(X_test_sc)

    # Compute ROC curve and ROC area
    fpr_mlp, tpr_mlp, _ = roc_curve(y_test, y_score_dnn[:,1])
    roc_auc_mlp = auc(fpr_mlp, tpr_mlp)
    
    print('ROC: Attempt {0} has an AUC of {1:.2f}.'.format(episode_count, roc_auc_mlp))
#    plt.figure(figsize=(13,8))
#    
#    plt.plot(fpr_mlp, tpr_mlp,
#         lw=2, label="ROC curve (AUC = {:.2f})".format(roc_auc_mlp))
#    
#    plt.rc('text', usetex=True)
#    plt.rc('font', family='serif')
#    plt.plot([0, 1], [0, 1], color='black', lw=2, linestyle='--')
#    plt.xlim([0.0, 1.0])
#    plt.ylim([0.0, 1.05])
#    plt.grid()
#    plt.xlabel('False Positive Rate')
#    plt.ylabel('True Positive Rate')
#    plt.title('Receiver operating characteristic -- Attempt {}'.format(episode_count))
#    plt.legend(loc="lower right")
#    plt.savefig('dnn_roc_{}.pdf'.format(episode_count), format='pdf')
    episode_count = episode_count + 1
    #plt.show()

    return [fpr_mlp, tpr_mlp, roc_auc_mlp] # model is valid
        
# create model
def create_mlp(width, depth, act):
    mlp = Sequential()
    mlp.add(Dense(units=width, input_dim=nX, activation=act))
    for k in np.arange(depth): # the depth of the neural network
        mlp.add(Dense(width, use_bias=False))

    mlp.add(Dense(units=nY, input_dim=width, activation=act))

    mlp.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])  # log loss is an alternative for AUC; ignore accuracy
    return mlp

#def log_loss(y_true, y_pred):
 #  return K.mean(K.binary_crossentropy(y_true, y_pred), axis=-1)