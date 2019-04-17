#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 12:22:25 2019

@author: Or Duek
Analysis with FSL
"""


from __future__ import print_function
from builtins import str
from builtins import range

import os  # system functions

from nipype import config
# config.enable_provenance()

from nipype.interfaces import spm, fsl

import scipy.io as spio
import nipype.interfaces.io as nio  # Data i/o
import nipype.interfaces.utility as util  # utility
import nipype.pipeline.engine as pe  # pypeline engine
import nipype.algorithms.rapidart as ra  # artifact detection
import nipype.algorithms.modelgen as model  # model specification

from nipype.algorithms.misc import Gunzip
from nipype import Node, Workflow, MapNode
from nipype import SelectFiles
from os.path import join as opj


os.chdir('/media/Data/work')
import readConditionFiles_r_aPTSD
# Specify the location of the data.
data_dir = os.path.abspath('/media/Data/FromHPC/output/fmriprep')

# Specify the subject directories
subject_list = ['1063' ]#, '1072'] # '1206', '1244','1273' ,'1291', '1305', '1340', '1345', '1346']
# Map field names to individual subject runs.
task_list = ['3','4','5','6']
infosource = Node(util.IdentityInterface(fields=['subject_id'
                                            ],
                                    ),
                  name="infosource")
infosource.iterables = [('subject_id', subject_list),('task_id',task_list)]

# Map field names to individual subject runs.
info = dict(
    func = [['subject_id', ['3', '4', '5', '6'],'preproc_bold']],
    mask= [['subject_id', ['3', '4', '5', '6'],'brain_mask']]
    )



infosource = Node(util.IdentityInterface(fields=['subject_id'
                                            ],
                                    ),
                  name="infosource")
infosource.iterables = [('subject_id', subject_list)]

# SelectFiles - to grab the data (alternativ to DataGrabber)
templates = {'func': '/media/Data/FromHPC/output/fmriprep/sub-%s/ses-1/func/sub-*_ses-1_task-%s_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz',
             'mask': '/media/Data/FromHPC/output/fmriprep/sub-%s/ses-1/func/sub-*_ses-1_task-%s_space-MNI152NLin2009cAsym_desc-brain_mask.nii.gz'}
#selectfiles = Node(SelectFiles(templates,
#                               base_directory=data_dir),
#                   name="selectfiles")

datasource = pe.Node(
    interface=nio.DataGrabber(
        infields=['subject_id'], outfields=['func','mask']),
    name='datasource')
datasource.inputs.base_directory = data_dir
datasource.inputs.template = 'sub-%s/ses-1/func/sub-*_ses-1_task-%s_space-MNI152NLin2009cAsym_desc-%s.nii.gz'
datasource.inputs.template_args = info
datasource.inputs.sort_filelist = True
#datasource.inputs.field_template = templates
os.chdir('/media/Data/work')

fsl.FSLCommand.set_default_output_type('NIFTI_GZ')

###########################
def subjectinfo(subject_id):
   
    from readConditionFiles_r_aPTSD import loadmat, _check_keys, _todict, readConditions
    import scipy.io as spio
    import os,json,glob,sys
    ###############################################################################################
    # Define experiment things (data_dir = filder where data files are present. )
    data_dir ='/media/Data/FromHPC/output/fmriprep'
    from bids.grabbids import BIDSLayout
    layout = BIDSLayout(data_dir)
    tasks = ['3','4','5','6'] # number of task (i.e. block. corresponding to file name)
    source_epi = layout.get(type="bold", session="1", extensions="nii.gz", subject = subject_id)
    ###############################################################################################
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
    
    
    # creates full table of subject info (conditions, runs etc.)
   # onsets = ([])
    import pandas as pd
    eventsTotal = organizeBlocks(subject_id)
    for i in range(len(eventsTotal)):
        print (i)
        eventsTotal[i]['condName'] = 'test'
        for n in range(1,len(eventsTotal[i])+1):
            if eventsTotal[i].condition[n] =='Gain':
                if eventsTotal[i].trial_type[n] == 'risk':
                    eventsTotal[i]['condName'][n] = 'GainRisk'
                else:
                    eventsTotal[i]['condName'][n] = 'GainAmb'
            if eventsTotal[i].condition[n] == 'Loss':
                if eventsTotal[i].trial_type[n] == 'risk':
                    eventsTotal[i]['condName'][n] = 'LossRisk'
                else:
                    eventsTotal[i]['condName'][n] = 'LossAmb'
    #events = eventsTotal[0]# the first ses
    
    from nipype.interfaces.base import Bunch
    #from copy import deepcopy
    print("Subject ID: %s\n" % str(subject_id))
    output = []
    #names = ['GainRisk', 'GainAmb','LossRisk', 'LossAmb']
    for r in range(4):
        
        print (r)
        confounds = pd.read_csv(os.path.join(data_dir, 
                                        "sub-%s"%subject_id, "ses-%s"%source_epi[r].session, "func", 
                                        "sub-%s_ses-%s_task-%s_desc-confounds_regressors.tsv"%(source_epi[r].subject, source_epi[r].session, tasks[r])),
                                           sep="\t", na_values="n/a")
        # put here ifs to build bunch for each run according to the conditions. 
        if eventsTotal[r].condition[r+1]=='Loss':
            print ('LOSS')
           
            output.insert(r,
                          Bunch(conditions = ['LossRisk','LossAmb','GainRisk' ,'GainAmb'],
                                onsets = [list(eventsTotal[r][eventsTotal[r].condName=='LossRisk'].onset),
                                         list(eventsTotal[r][eventsTotal[r].condName=='LossAmb'].onset),
                                         [0],[0]],
                                durations = [list(eventsTotal[r][eventsTotal[r].condName=='LossRisk'].duration),
                                             list(eventsTotal[r][eventsTotal[r].condName=='LossAmb'].duration),
                                             [0],[0]
                                             ],
                                             regressors=[list(confounds.framewise_displacement.fillna(0)),
                                                         list(confounds.a_comp_cor_00),
                                                         list(confounds.a_comp_cor_01),
                                                         list(confounds.a_comp_cor_02),
                                                         list(confounds.a_comp_cor_03),
                                                         list(confounds.a_comp_cor_04),
                                                         list(confounds.a_comp_cor_05),
                                                         ],
                             regressor_names=['FramewiseDisplacement',
                                              'aCompCor0',
                                              'aCompCor1',
                                              'aCompCor2',
                                              'aCompCor3',
                                              'aCompCor4',
                                              'aCompCor5'],
                         
                      
                               
                                                                  ) )
        elif eventsTotal[r].condition[r+1]=='Gain':
            print ('Gain')
            
            output.insert(r,
                          Bunch(conditions = ['LossRisk','LossAmb','GainRisk','GainAmb'],
                                onsets = [[0],[0],
                                         list(eventsTotal[r][eventsTotal[r].condName=='GainRisk'].onset),
                                         list(eventsTotal[r][eventsTotal[r].condName=='GainAmb'].onset)],
                                durations = [[0],[0],
                                             list(eventsTotal[r][eventsTotal[r].condName=='GainRisk'].duration),
                                             list(eventsTotal[r][eventsTotal[r].condName=='GainAmb'].duration)],
                                             regressors=[list(confounds.framewise_displacement.fillna(0)),
                         list(confounds.a_comp_cor_00),
                         list(confounds.a_comp_cor_01),
                         list(confounds.a_comp_cor_02),
                         list(confounds.a_comp_cor_03),
                         list(confounds.a_comp_cor_04),
                         list(confounds.a_comp_cor_05),
                        ],
             regressor_names=['FramewiseDisplacement',
                              'aCompCor0',
                              'aCompCor1',
                              'aCompCor2',
                              'aCompCor3',
                              'aCompCor4',
                              'aCompCor5'],
                        
             
                              )
                               
                            
                                )
    return output


#########################################
#Adding subjectInfo function as a node

# Get Subject Info - get subject specific condition information
getsubjectinfo = Node(util.Function(input_names=['subject_id'],
                               output_names=['subject_info'],
                               function=subjectinfo),
                        
                      name='getsubjectinfo')

################################################################

# creating contrasts
condition_names = ['GainRisk', 'GainAmb' ,'LossRisk', 'LossAmb']                          

GainRisk_cond = ['GainRisk','T', condition_names ,[1,0,0,0]]
GainAmb_cond = ['GainAmb','T', condition_names ,[0,1,0,0]]
LossRisk_cond = ['LossRisk','T', condition_names,[0,0,1,0]]
LossAmb_cond = ['LossAmb','T',condition_names,[0,0,0,1]]
Risk_vs_Amb = ["Risk vs. Amb",'T', condition_names ,[0.5, -0.5, 0.5, -0.5]]
Gain_vs_Loss = ["Gain vs. Loss",'T', condition_names ,[0.5, 0.5, -0.5, -0.5]]

#all_motor = ["All motor", 'F', [LossRisk_cond, LossAmb_cond, GainRisk_cond, GainAmb_cond]]
contrasts=[GainRisk_cond, GainAmb_cond, LossRisk_cond, LossAmb_cond, Risk_vs_Amb, Gain_vs_Loss]
#%% merge all images to one file
merge = Node(interface = fsl.Merge(), name = 'merge')
merge.inputs.dimension = 't'
testWf = Workflow(name="testingMergeFsk", base_dir="/media/Data/work")                   
testWf.connect([
        (infosource, datasource, [('subject_id','subject_id')]),
        (datasource, merge, [('func','in_files')]),
        ])

testWf.run()
#%%
os.chdir('/media/Data/work')
skip = MapNode(interface = fsl.ExtractROI(), name = "modelROIext", base_dir = '/media/Data/work', iterfield=['in_file'])
skip.inputs.t_min = 4
skip.inputs.t_size = -1
skip.inputs.output_type = 'NIFTI_GZ'
#skip.inputs.in_file = '/media/Data/FromHPC/output/fmriprep/sub-1063/ses-1/func/sub-1063_ses-1_task-3_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz'
#h= skip.run()

modelspec = Node(interface=model.SpecifyModel(), name="modelspec") 

modelspec.inputs.input_units = 'scans'
#modelspec.inputs.output_units = 'scans'
#modelspec.inputs.outlier_files = '/media/Data/R_A_PTSD/preproccess_data/sub-1063_ses-01_task-3_bold_outliers.txt'
modelspec.inputs.time_repetition = 1  # make sure its with a dot 
modelspec.inputs.high_pass_filter_cutoff = 128.



level1design = Node(interface=fsl.model.Level1Design(), name = "level1Design")
level1design.inputs.interscan_interval = 1
level1design.inputs.bases = {'dgamma':{'derivs': True}}
level1design.inputs.model_serial_correlations=True
level1design.inputs.contrasts=contrasts

modelgen = MapNode(interface=fsl.model.FEATModel(), name = "modelgen", iterfield = ['fsf_file','ev_files'])
# in_file is fsf_file=level1design_results.outputs.fsf_files,  ev_files=level1design_results.outputs.ev_file

mask = MapNode(interface=fsl.maths.ApplyMask(), name = "mask", iterfield =['mask_file', 'in_file'])
# in_file is subject info

filmgls= MapNode(interface=fsl.FILMGLS(), name = "filmgls", iterfield=['design_file', 'in_file', 'tcon_file'])
filmgls.inputs.autocorr_noestimate = True
filmgls.inputs.results_dir = 'stats'



'''
in_file=mask_results.outputs.out_file,
design_file = modelgen_results.outputs.design_file,
tcon_file = modelgen_results.outputs.con_file,
fcon_file = modelgen_results.outputs.fcon_file,
'''

wFSL = Workflow(name="l1FSL", base_dir="/media/Data/work")                       



wFSL.connect([
        (infosource, datasource, [('subject_id','subject_id')]),
        (datasource, skip, [('func','in_file')]),
        (skip, modelspec, [('roi_file', 'functional_runs')]),
        (infosource, getsubjectinfo, [('subject_id', 'subject_id')]),
        (getsubjectinfo,modelspec, [('subject_info', 'subject_info')]),
        (skip, mask, [('roi_file','in_file')]),
        (datasource, mask, [('mask','mask_file')])
        
        ])
wFSL.connect([(modelspec, level1design, [("session_info", "session_info")])])


wFSL.connect([
         (level1design, modelgen, [('fsf_files','fsf_file'), ('ev_files', 'ev_files')]),
         (mask, filmgls, [('out_file', 'in_file')]),
         (modelgen, filmgls, [('design_file','design_file'),('con_file','tcon_file'),('fcon_file','fcon_file')]),
         
    ])

#######
#%% merging map of filmgls
###################
copemerge = pe.MapNode(
    interface=fsl.Merge(dimension='t'),
    iterfield=['in_files'],
    name="copemerge")

varcopemerge = pe.MapNode(
    interface=fsl.Merge(dimension='t'),
    iterfield=['in_files'],
    name="varcopemerge")                   

#%%
wFSL.connect([
        (filmgls, copemerge, [('copes','in_files')]),
        (filmgls, varcopemerge, [('varcopes','in_files')])
        ])
#%% 2nd level
level2model = pe.Node(interface=fsl.L2Model(), name='l2model')
level2model.inputs.num_copes = 6

flameo = pe.MapNode(
    interface=fsl.FLAMEO(run_mode='fe'),
    name="flameo",
    iterfield=['cope_file', 'var_cope_file'])

#%%
def sort_copes(files):
    numelements = len(files[0])
    outfiles = []
    for i in range(numelements):
        outfiles.insert(i, [])
        for j, elements in enumerate(files):
            outfiles[i].append(elements[i])
    return outfiles


def num_copes(files):
    return len(files)


#%% adding second level to workflow
wFSL.connect([
    (copemerge, flameo, [('merged_file', 'cope_file')]),
    (varcopemerge, flameo, [('merged_file', 'var_cope_file')]),
    (level2model, flameo, [('design_mat', 'design_file'),
                           ('design_con', 't_con_file'), ('design_grp',
                                                          'cov_split_file')]),
])


#%%
wFSL.write_graph("workflow_graph.dot", graph2use='colored', format='png', simple_form=True)
from IPython.display import Image
Image(filename="/media/Data/work/l1FSL/workflow_graph.png")
wFSL.write_graph(graph2use='flat')

Image(filename = '/media/Data/work/l1run/graph_detailed.png')

wFSL.run('MultiProc', plugin_args={'n_procs': 3})         


#% Graph single subject
from nilearn.plotting import plot_stat_map
from nilearn.plotting import plot_glass_brain
anatimg = '/media/Data/FromHPC/output/fmriprep/sub-1072/anat/sub-1072_space-MNI152NLin2009cAsym_desc-preproc_T1w.nii.gz'
plot_stat_map(
    '/media/Data/work/l1FSL/_subject_id_1072/filmgls/mapflow/_filmgls1/stats/zstat1.nii.gz', title='LossRisk - fwhm=6',
    bg_img=anatimg, threshold=3, display_mode='x', cut_coords=(-5, 0, 5, 10, 15), dim=0);
        

plot_glass_brain(
    '/media/Data/work/l1FSL/_subject_id_1072/filmgls/mapflow/_filmgls0/stats/cope1.nii.gz', colorbar=True,
     threshold=3.6, display_mode='lyrz', black_bg=True, vmax=10, title='spm_fwhm6_ GainRIsk vs Baseline');       

plot_stat_map(
    '/media/Data/work/l1FSL/_subject_id_1072/filmgls/mapflow/_filmgls0/stats/tstat3.nii.gz',
     threshold=3.6, display_mode='ortho');