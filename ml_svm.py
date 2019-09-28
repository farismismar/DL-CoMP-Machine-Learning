#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 28 17:48:29 2019

@author: farismismar
"""

import os
###################################################

#os.chdir('/Users/farismismar/Desktop/DL CoMP/')


from sklearn.svm import SVC
from sklearn.model_selection import GridSearchCV

import pandas as pd
from sklearn.model_selection import train_test_split
import numpy as np
import random

from sklearn.metrics import roc_auc_score, accuracy_score

from scipy.io import loadmat

# Set the random seed
seed = 0

# save X and Y dimensions   
mX, nX = [0,0]
mY, nY = [0,0]

y_pred = None # placeholder
model = None # the forward definition.

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

    X_test.fillna(0, inplace=True, axis=0)
    #print(X_test)
    try:
        y_pred = model.predict(X_test)
        print('INFO: Prediction is y = {}'.format(y_pred[0]));
    except:
        #print('WARNING: Invalid prediction using {}'.format(X_test.values));
        y_pred = np.nan

    return y_pred

# returns the test error
def train_wrapper(filename): # filename='measurements.mat'
    global dataset
    global seed
    global model
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

    # Perform a split 30-70
    train, test = train_test_split(dataset, test_size=0.30, random_state=seed)
    
    X_train = train.drop('y', axis=1)
    X_test = test.drop('y', axis=1)
    
    y_train = train['y'].values
    y_test = test['y'].values

    mX, nX = X_train.shape
    mY = y_train.shape
    nY = 1

    # The hyperparameters
    p = 4
    Cs = [0.01, 1, 10]
    kernels = ['rbf', 'linear', 'poly']
    param_grid = {'C': Cs, 'kernel': kernels}
    grid_search = GridSearchCV(SVC(degree=p, gamma='scale', class_weight='balanced', probability=True), param_grid, verbose=True, n_jobs=1, cv=3) # p-degree kernel.

    try:
        grid_result = grid_search.fit(X_train, y_train)
    except:
        print('WARNING: Single class classifier.  Fallback.')
        return [np.nan, np.nan]

    # This is the best model
    best_model_svm = grid_result.best_params_
    print(best_model_svm)
 
    model = grid_result.best_estimator_
    svc = model

    y_pred = svc.predict(X_test)
    y_score = svc.predict_proba(X_test)

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
        

# Local code for debugging
#initialize_wrapper(seed)
#train_wrapper('measurements.mat')
