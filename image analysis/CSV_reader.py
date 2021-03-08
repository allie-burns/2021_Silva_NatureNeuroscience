#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  9 15:45:15 2020

@author: lukasvandenheuvel

This script calculates 
(1) the ratio of cFos positive cells, 
(2) the ratio of cFos and tracer positive (double positive) cells, and 
(3) the ratio of double-positive cells normalized by chance, 
for specified brain regions in analyzed microscopy slices.

INPUT
-----
A set of .csv files with cell counts (cFos+, tracer+ and dapi+) on microscopy images.
These files are generated using the QuPath software, with help of the cell-count script 
(qps-multiChannels_PositiveCells_Analysis.groovy). There is one file for each microscopy image.

OUTPUTS
-------
(1) A .csv file called "experiment + '_results_' + region_name + '_' + tracer + '.csv'" 
inside the specified output directory. Example: WT47_results_INS_BLA.csv
This file contains the ratios (cFos, double-positive, and double-positive / chance),
for each animal and each specified brain region, averaged over the brain slices available for each animal.
The standard errors of the mean (SEM) are also provided.

(2) A folder called "Raw_numbers'+'_'+region_name+'_'+tracer" containing one .csv file per animal. 
These files contain the cell counts used to calculate the three ratios.

This code was used to perform the analysis presented in "A thalamo-amygdalar circuit underlying the extinction 
of remote fear memories" (Silva et al., Nature Neuroscience, 2021).
"""

# import functions in the file QuPath.py:
from QuPath import import_files, calculate_ratio, get_chance_ratio, write_results_to_csv, output_raw_numbers

#%%
'''----------------------------PARAMETERS THAT YOU HAVE TO SET----------------------------'''

# specifications
experiment = 'WT48'      # Experiment name.
region_name = 'mPFC'     # Is only used in the title of the output file. Choose whatever you like.
tracer = 'BLA'           # 'BLA' or 'NRe'

# where are the .csv files stored?
root = './example data/cell counts in mPFC'

# Path and name of the directory where you want to output. If it is a path to a folder that doesn't exist yet, it will be created automatically.
output_directory = './example data/ratios in mPFC' 

# which mice are in which groups?
recall = []
# WT48
extinction = ['14031', '14032', '14033', '14034', '14035', '14036']
control = ['14025', '14026', '14027', '14028', '14029', '14030']

# WT47
#extinction = ['12852', '12853', '12854', '12855']
#control = ['12854(50)', '12851', '12857', '12858']

# What are the whole regions and subregions? Keys are whole regions, values are lists of subregions.
# Example:

'''
# Midline thalamus
combine = {'NRe': ['Re', 'Xi', 'PaXi'],
           'aNRe': ['aRe', 'aXi', 'aPaXi'],
           'VRe': ['VRe', 'aVRe'],
           'PVT': ['PVT'],
           'IMD': ['IMD'],
           'MD': ['MD'],
           'CM': ['CM'],
           'aCM': ['aCM'],
           'IAM': ['IAM'],
           'Rh': ['Rh'],
           'aRh': ['aRh'],
           'Patch': ['aPatch'],
           'pPatch': ['pPatch']
          }
'''

'''
combine = {'Re': ['Re', 'aRe'],
           'VRe': ['VRe', 'aVRe'],
           'Xi': ['Xi', 'aXi'],
           'PaXi': ['PaXi', 'aPaXi'],
           'PVT': ['PVT'],
           'IMD': ['IMD'],
           'MD': ['MD'],
           'CM': ['CM', 'aCM'],
           'IAM': ['IAM'],
           'Rh': ['Rh', 'aRh'],
           'Patch': ['aPatch', 'pPatch']
          }
'''
'''
combine = {'Re': ['Re'],
           'aRe': ['aRe'],
           'VRe': ['VRe'],
           'aVRe': ['aVRe'],
           'Xi': ['Xi'],
           'aXi': ['aXi'],
           'PaXi': ['PaXi'], 
           'aPaXi': ['aPaXi'],
           'PVT': ['PVT'],
           'IMD': ['IMD'],
           'MD': ['MD'],
           'CM': ['CM'],
           'aCM': ['aCM'],
           'IAM': ['IAM'],
           'Rh': ['Rh'], 
           'aRh': ['aRh'],
           'pPatch': ['pPatch'],
           'aPatch': ['aPatch'] 
          }
'''

'''

combine = {'NRe': ['Re', 'Xi', 'PaXi', 'aRe', 'aXi', 'aPaXi'],
           'VRe': ['VRe', 'aVRe'],
           'PVT': ['PVT'],
           'IMD': ['IMD'],
           'MD': ['MD'],
           'CM': ['CM', 'aCM'],
           'IAM': ['IAM'],
           'Rh': ['Rh', 'aRh'],
           'Patch': ['aPatch', 'pPatch'],
          }
'''

'''
# VTA
combine = {'VTA': ['VTA'],
           'mVTA': ['mVTA'],
           'SUM': ['SUM'],
           'IF': ['IF'],
           'LN': ['LN'],
           'ZI': ['ZI'],
           'IPN': ['IPN'],
           'PBP': ['PBP'],
           'IPR': ['IPR'],
           'IPL': ['IPL'],
           'IPI': ['IPI'],
           'IPC': ['IPC'],
           'IPF': ['IPF'],
           'PIF': ['PIF']
        }
'''

# mPFC
combine = {'IL': ['IL23', 'IL5', 'IL6'],
           'PL': ['PL23', 'PL5', 'PL6'],
           'AC': ['AC23', 'AC5', 'AC6'],
          }
'''
# Insular
combine = {'GI': ['GI'],
           'DI': ['DI'],
           'AI': ['AI']
          }
'''

#%%
'''-------------------------------START CODE------------------------------------'''

# Output files
output_file = experiment + '_results_' + region_name + '_' + tracer + '.csv'    # name of the output file
title = experiment + ' ' + region_name + '->' + tracer  # title (will appear on the first line of the output file)

# Lists of parameters, whole regions and all mice
param_list = ['DAPI','cFos',tracer,'area',tracer + '_cFos']
whole_regs = list(combine.keys())
full_mouse_list = recall + extinction + control 

# only consider regions where we have at least x amount of traced cells:
traced_cell_thresh_VMT = 5
traced_cell_thresh_single_region = 4

#%% import all files in the channel folders into the mouse dictionary and, if necessary, combine regions.
m_dict = import_files(root,combine,param_list)
 
#%% analyse cFos/tracer cell counts
[cFos_ratio_df, cFos_SEM_df] = calculate_ratio('cFos', 'DAPI', m_dict, full_mouse_list, whole_regs, traced_cell_thresh_single_region)
[Doublepos_ratio_df, Doublepos_SEM_df] = calculate_ratio(tracer+'_cFos', tracer, m_dict, full_mouse_list, whole_regs, traced_cell_thresh_single_region)

#%% tracer+cFos chance ratio
[BLA_chance_df,BLA_chance_SEM_df] = get_chance_ratio(tracer, m_dict, full_mouse_list, whole_regs, traced_cell_thresh_single_region) 

#%% Write output to csv file
write_results_to_csv(output_directory, output_file, title, recall, extinction, control, tracer,
                     BLA_chance_df, BLA_chance_SEM_df, 
                     cFos_ratio_df, cFos_SEM_df, 
                     Doublepos_ratio_df, Doublepos_SEM_df)

output_raw_numbers(output_directory, m_dict, tracer, region_name)
