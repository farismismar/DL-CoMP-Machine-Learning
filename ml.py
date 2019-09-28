# -*- coding: utf-8 -*-
"""
Created on Thu April 28 20:46:45 2017

@author: farismismar
"""
# Use these lines to force GPU
import os
os.environ["CUDA_DEVICE_ORDER"]="PCI_BUS_ID"
# The GPU id to use, usually either "0" or "1"
os.environ["CUDA_VISIBLE_DEVICES"]="0" 
###################################################

#os.chdir('/Users/farismismar/Desktop/DL CoMP/')

# Do other imports now...
import tensorflow as tf

from keras.models import Sequential
from keras.layers import Dense
from keras.optimizers import SGD

from sklearn.utils import class_weight
from sklearn.model_selection import GridSearchCV, train_test_split
from keras.wrappers.scikit_learn import KerasClassifier
from keras.backend.tensorflow_backend import set_session

import pandas as pd
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
#from matplotlib import rc
import numpy as np
import random

from sklearn.preprocessing import MinMaxScaler

from sklearn.metrics import roc_auc_score, accuracy_score

from scipy.io import loadmat

# Set the random seed
seed = 0
tf.set_random_seed(seed)
np.random.seed(seed)
random.seed(seed)

# save X and Y dimensions   
mX, nX = [0,0]
mY, nY = [0,0]

y_pred = None # placeholder

model = None # the forward definition.
ss = None # The scaler forward def

def initialize_wrapper(random_state):
    global seed
    seed = random_state
    np.random.seed(random_state)
    random.seed(random_state)

def predict_wrapper(RSRP, SINR_1, SINR_2, rank):
    X_test = pd.DataFrame({'RSRP': RSRP,
                           'TBSINR_1': SINR_1,
                           #'TBSINR_2': SINR_2,
                            'rank' :rank}, index=[0])

    global model
    global y_pred
    global ss

    X_test.fillna(0, inplace=True, axis=0)

    try:
        y_pred = model.predict(ss.transform(X_test))
        print('INFO: Prediction is y = {}'.format(y_pred[0]));
    except:
        #print('WARNING: Invalid prediction.');
        y_pred = np.nan

    return y_pred

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
     
    data = loadmat(filename) # this is a dict.
    keys = list(data.keys())[3:] # skip the first three columns
    values = list(data.values())[3:]
    
    dataset = pd.DataFrame()
    dataset = dataset.reindex(columns = keys) # create an empty dataframe

    for ii in np.arange(len(values)):
        v_ = np.array(values[ii])
        dataset[keys[ii]] = pd.Series(v_.flatten()) # cannot add the data to this empty df.
    
    dataset['y'] = 1*(dataset['BLER'] <= 0.1) # H-ARQ target.
    dataset = dataset[['RSRP', 'TBSINR_1', 'rank', 'y']]
    dataset.dropna(inplace=True, axis=0)
    if os.path.exists('dataset.csv'):
        dataset.to_csv('dataset.csv', index=False, mode='a', header=False) # append
    else:
        dataset.to_csv('dataset.csv', index=False)

    #print(dataset.head())
   
    # Perform a split 30-70
    train, test = train_test_split(dataset, test_size=0.30, random_state=seed)
    
    X_train = train.drop('y', axis=1)
    X_test = test.drop('y', axis=1)
    
    y_train = train['y'].values
    y_test = test['y'].values

    mX, nX = X_train.shape
    mY = y_train.shape
    nY = 1
    
    ss = MinMaxScaler(feature_range=(0,1))
    
    # Scale the variables
    X_train_sc = ss.fit_transform(X_train)
    X_test_sc = ss.transform(X_test)
    
    model = KerasClassifier(build_fn=create_mlp, verbose=0, epochs=10, batch_size=8)

    # The hyperparameters
    width_dims=[3,5,10]
    n_hiddens = [3,5] # the depth of hidden layers

    hyperparameters = dict(width=width_dims, depth=n_hiddens) 
    class_weights = class_weight.compute_class_weight('balanced', np.unique(y_train), y_train)

    grid = GridSearchCV(estimator=model, param_grid=hyperparameters, n_jobs=1, cv=3)
    gpu_available = tf.test.is_gpu_available()
    if (gpu_available == False):
        print('WARNING: No GPU available.  Will continue with CPU.')
        
    with tf.device('/gpu:0'):
        grid_result = grid.fit(X_train_sc, y_train, class_weight=class_weights)
    
    # This is the best model
    best_model_mlp = grid_result.best_params_
    print(best_model_mlp)
 
    model = grid_result.best_estimator_
    mlp = model

    y_pred = mlp.predict(X_test_sc)
    y_score = mlp.predict_proba(X_test_sc)

    mu = accuracy_score(y_test, y_pred)
    # Compute ROC curve and ROC area
    try:
        roc_auc = roc_auc_score(y_test, y_score[:,1])    
    except:
        print('WARNING: ROC was not computed.  Returning NaN');
        roc_auc = np.nan
 
    print('ROC for training is: {}'.format(roc_auc))
    print('Misclassification error for training is: {:.3f}'.format(1-mu))
    
    return [roc_auc, 1-mu] # model is valid
        
# create model
def create_mlp(width, depth):
    mlp = Sequential()
    mlp.add(Dense(units=width, input_dim=nX, use_bias=True, activation='sigmoid'))
    for k in np.arange(depth): # the depth of the neural network
        mlp.add(Dense(width))
    mlp.add(Dense(units=nY, input_dim=width, use_bias=True, activation='sigmoid'))

    mlp.compile(loss='binary_crossentropy', optimizer=SGD(lr=0.01), metrics=['accuracy'])  # log loss is an alternative for AUC; ignore accuracy
    return mlp

# Local code for debugging
#initialize_wrapper(seed)
#train_wrapper('measurements.mat')
