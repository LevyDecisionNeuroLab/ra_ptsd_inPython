# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 00:06:31 2019

@author: ruonanjia
"""

import scipy.io as spio
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# change working directory
os.chdir('D:\\Ruonan\\Projects in the lab\\VA_RA_PTB\\Imaging analysis\\imaging_analysis_inpython_041719\\ra_ptsd_inPython\\')
# load functions
from readConditionFiles_r_aPTSD import loadmat, _check_keys, _todict


#%% functions borrowed from Or and modified
def organizeBlocks(subNum, path):
    # Read both mat files (first timestamp)
    # check first block of each day. 
    # check thrird of each day
    # sort
    orderArray = []
    matFileLoss = os.path.join(path, 'RA_LOSS_%s_fitpar.mat'%subNum) 
    matFileGain = os.path.join(path, 'RA_GAINS_%s_fitpar.mat'%subNum) 
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

#%%
def readConditions(matFile, x): # takes name of file and when to begin (i.e. first block is zero. second is 32 etc.)
    metaData = loadmat(matFile)
    timeStamp = []
    condition = []
    events = []
    resultsArray = []
    duration = []
    sv = []
    sv_label_val = []
    sv_label_sal = []
    
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
       
       sv.append(metaData['Data']['sv'][i])
       sv_label_val.append(metaData['Data']['sv_label_val'][i])
       sv_label_sal.append(metaData['Data']['sv_label_sal'][i])
    
    events= pd.DataFrame({'trial_type':condition, 'onset':timeStamp, 'duration':duration, 'sv':sv, 'sv_label_val':sv_label_val, 'sv_label_sal':sv_label_sal})[1:] # building data frame from what we took. Removing first row because its not used. 
    return events
    

#%%
root = 'D:\\Ruonan\\Projects in the lab\\VA_RA_PTB\\'

fitpar_dir = os.path.join(root, 'Analysis Ruonan\\Fitpar files\\Behavior data fitpar_04232019')


mat_filename = os.path.join(fitpar_dir, 'RA_GAINS_3_fitpar.mat')

sub_data = loadmat(mat_filename)
sub_data.keys()


data=sub_data['Data']    
for key in data:
    print(key)
    
trial_time_meta = data['trialTime'] 

trial_time = vars(trial_time[0]) # a dict for trial 0

trial_time['trialStartTime'] # array (6,): year,month,day,hr,min,sec

ambigs = data['ambigs']

sv_sal = data['sv_label_sal']
# find trial start time for every trial, plus 6tr, and take them out
# label them based on sv-high/low, domain-gain/loss

#%% create mask and labels for a subject, all functional runs

# read events for a single run
events = readConditions(mat_filename,)
# read events for all (8) runs
total_event = organizeBlocks(1069, fitpar_dir)

# products: label of the whole 500 trs, mask of trs to take out the 6 TR duration
tr_mask = [] # each list element is a run, length should equals to 500
sv_val_label = []
sv_sal_label = []
sv = []

for event_run in total_event:
    sv_run = np.zeros(500)
    tr_mask_run = np.zeros(500)
    sv_val_label_run = np.zeros(500)
    sv_sal_label_run = np.zeros(500)
    
    # onsets in TR
    onset = event_run.onset
    
    # create mask
    for tr_delta in range(0,event_run.loc[1, 'duration']):
        tr_mask_run[onset+tr_delta] = 1
        sv_run[onset+tr_delta] = event_run.loc[:,'sv']
        sv_val_label_run[onset+tr_delta] = event_run.loc[:, 'sv_label_val']
        sv_sal_label_run[onset+tr_delta] = event_run.loc[:, 'sv_label_sal']
    
    tr_mask.append(tr_mask_run)    
    sv.append(sv_run)
    sv_val_label.append(sv_val_label_run)
    sv_sal_label.append(sv_sal_label_run)


# check by plotting    
plt.plot(tr_mask[0])

plt.plot(sv_val_label[0][tr_mask[0]==1])
plt.plot(sv_sal_label[0][tr_mask[0]==1])     
plt.plot(sv[0][tr_mask[0]==1])  
                   
        
                   
    
    
    
    


