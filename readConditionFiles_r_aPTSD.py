#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 10:48:03 2019

@author: or

Loading and reading mat file
"""

import scipy.io as spio

def loadmat(filename):
    '''
    this function should be called instead of direct spio.loadmat
    as it cures the problem of not properly recovering python dictionaries
    from mat files. It calls the function check keys to cure all entries
    which are still mat-objects
    
    from: `StackOverflow <http://stackoverflow.com/questions/7008608/scipy-io-loadmat-nested-structures-i-e-dictionaries>`_
    '''
    data = spio.loadmat(filename, struct_as_record=False, squeeze_me=True)
    return _check_keys(data)

def _check_keys(dict):
    '''
    checks if entries in dictionary are mat-objects. If yes
    todict is called to change them to nested dictionaries
    '''
    for key in dict:
        if isinstance(dict[key], spio.matlab.mio5_params.mat_struct):
            dict[key] = _todict(dict[key])
    return dict        

def _todict(matobj):
    '''
    A recursive function which constructs from matobjects nested dictionaries
    '''
    dict = {}
    for strg in matobj._fieldnames:
        elem = matobj.__dict__[strg]
        if isinstance(elem, spio.matlab.mio5_params.mat_struct):
            dict[strg] = _todict(elem)
        else:
            dict[strg] = elem
    return dict

import pandas as pd

def readConditions(matFile, x): # takes name of file and when to begin (i.e. first block is zero. second is 32 etc.)
    metaData = loadmat(matFile)
    timeStamp = []
    condition = []
    events = []
    resultsArray = []
    duration = []
    #x= 0 # where to start (each block is 32)
    ambigs = metaData['Data']['ambigs']
    for i in range(x, x+31): # every block is 32 trials. Thus we stop at 31
        #a= metaData['Data']['trialTime'][i]
       # b = vars(a)
      
       resultsArray = vars(metaData['Data']['trialTime'][i])['trialStartTime'] - vars(metaData['Data']['trialTime'][x])['trialStartTime']
       timeStamp.append(int(round((3600*resultsArray[3] + 60*resultsArray[4] + resultsArray[5] -9)))) # using int and round to round to the close integer. 
       duration.append(6)
       
       if ambigs[i] == 0:
           condition.append('risk')
       else:
           condition.append('amb')
    
    events= pd.DataFrame({'trial_type':condition, 'onset':timeStamp, 'duration':duration})[1:] # building data frame from what we took. Removing first row because its not used. 
    return events



def organizeBlocks(subNum):
    # Read both mat files (first timestamp)
    # check first block of each day. 
    # check thrird of each day
    # sort
    orderArray = []
    matFileLoss = '/media/Drobo/Levy_Lab/Projects/R_A_PTSD_Imaging/Data/Behavior data/Behavior_fitpar/Behavior data fitpar_091318/RA_LOSS_%s_fitpar.mat'%subNum
    matFileGain = '/media/Drobo/Levy_Lab/Projects/R_A_PTSD_Imaging/Data/Behavior data/Behavior_fitpar/Behavior data fitpar_091318/RA_GAINS_%s_fitpar.mat'%subNum
    metaDataLoss = loadmat(matFileLoss)
    metaDataGain = loadmat(matFileGain)
    a= {'3rdLoss':list(vars(metaDataLoss['Data']['trialTime'][62])['trialStartTime']), '1stLoss':list(vars(metaDataLoss['Data']['trialTime'][0])['trialStartTime']), '1stGain':list(vars(metaDataGain['Data']['trialTime'][0])['trialStartTime']), '3rdGain':list(vars(metaDataGain['Data']['trialTime'][62])['trialStartTime'])}
    s = [(k, a[k]) for k in sorted(a, key=a.get, reverse=False)]
    for k, v in s:
         print (k, v)
         orderArray.append(k)
    totalEvent = []
    for n in orderArray:
        print (n)
        if n=='1stLoss':
            # run loss mat file on redConcitions function on first two blocks (i.e. 0, 31)
            print (n)
            for x in [0,31]:
                event = readConditions(matFileLoss, x)
                event['condition'] = 'Loss'
                totalEvent.append(event)
        elif n=='1stGain':
            # run Gain mat file on reCondition function
            print (n)
            for x in [0,31]:
                event = readConditions(matFileGain, x)
                event['condition'] = 'Gain'
                totalEvent.append(event)
        elif n=='3rdLoss':
            # run loss from 3rd block (i.e. 62, 93)
            print (n)
            for x in [62, 93]:
                event = readConditions(matFileLoss, x)
                event['condition'] = 'Loss'
                totalEvent.append(event)
        elif n=='3rdGain':
            # run gains from 3rd block
            print (n)
            for x in [62, 93]:
                event = readConditions(matFileGain, x)
                event['condition'] = 'Gain'
                totalEvent.append(event)
        else:
            print ('The condition ' + n + ' is not clear.')
        
        # the end result is an array of data sets per each run (i.e. block) - called totalEvent
    return totalEvent
    
    
    









