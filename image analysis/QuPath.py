#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 07:16:32 2020

@author: lukasvdh
"""

import numpy as np
import pandas as pd
import csv
import os
import shutil

#%%

def init_dict(mouse_list):
    '''
    This function initializes an empty dictionary with mouse names as keys.
    The dictionary will later be filled with dataframes (one dataframe per slice).
    ''' 
    m_dict = {}
    for mouse_name in mouse_list:
        m_dict[mouse_name] = []
    return m_dict

#%%

def find_regs_in_current_slice(df, combine):
    '''
    This function finds:
    - the regions that are present in the current slice;
    - which of the whole regions and subregions are present.
    The current regions are put in a list curr_regs 
    and whole regions with corresponding subregions are put in a dictionary curr_whole_regs.
    
    Example: 
    The regions present in the .csv file of the slice are 'Re', 'Xi' and 'MD left', 
    and Re and Xi must be combined into NRe (as specified in 'combine').
    Then curr_regs = ['Re', 'Xi', 'MD'],
    and curr_whole_regs =  {'NRe': ['Re', 'Xi']
                            'MD': ['MD']}
    '''
    curr_regs = []
    valid_regs = []
    curr_whole_regs = {}
    
    # find valid regions (the subregions that appear in 'combine')
    for whole_reg, subregs in combine.items():
        for subreg in subregs:
            valid_regs.append(subreg)
    
    # Find current regions:
    for reg in df.index:
        region = reg.split()[0]
        curr_regs.append(region)
        if reg.split()[0] not in valid_regs: # check if region in csv file is valid
            print('WARNING: csv file contains invalid region {r}!'.format(r = reg))
    # Remove doubles:
    curr_regs = list(dict.fromkeys(curr_regs))
    
    # Find current whole regions and subregions:
    for whole_reg, subregs in combine.items():
        if len(list(set(subregs) & set(curr_regs)))!= 0: 
            curr_whole_regs[whole_reg] = []
        for subreg in subregs:
            if subreg in curr_regs:
                curr_whole_regs[whole_reg].append(subreg)
    
    return curr_regs, curr_whole_regs

#%%

def find_current_mice(file_list): 
    '''
    This function makes a list of all mice that have a file in the dictionary 'root'.
    '''
    curr_mouse_list = []
    for file_name in file_list:
        curr_mouse_list.append(file_name.split('_')[0])

    curr_mouse_list = list(dict.fromkeys(curr_mouse_list))
    return curr_mouse_list

#%%
def csv_to_dataframe(data, param_list):
    '''
    This function takes as input raw data from a csv file (data = a dataframe created with pd.read_csv).
    It converts counts (Num A, Num B, etc) to number of detected cells ('DAPI, cFos', etc).
    To do this it sums the relevant counts, depending on what kind of analysis was done (1, 2 or 3 channel).
    
    param_list is a list with relevant parameters, e.g. ['DAPI, cFos, BLA, BLA_cFos']
    '''    
    # Make an empty table with rows = all regions in current slice, columns = param_list
    df = pd.DataFrame(np.nan, index=data.index, columns=param_list)
    
    # Determine which channels are counted:
    channels = data.columns 
        
    # Fill the dataframe
    #if 'Num C' not in channels and 'Num B' not in channels:  # 1 channel analysis (cFos only)
    if 'BLA' not in param_list and 'NRe' not in param_list:
        df['DAPI'] = data['Num Detections']
        df['cFos'] = data['Num A'] 
        df['area'] = data['Area um^2'] / 1e6
    #elif 'Num C' not in channels and 'Num B' in channels:    # 2 channel analysis (cFos+Re) with channel B --> NRe
    elif 'NRe' in param_list and 'BLA' not in param_list:
        df['DAPI'] = data['Num Detections']
        df['cFos'] = data['Num A'] + data['Num AB']
        df['NRe'] = data['Num B'] + data['Num AB']
        df['NRe_cFos'] = data['Num AB']
        df['area'] = data['Area um^2'] / 1e6
    #elif 'Num C' in channels and 'Num B' not in channels:    # 2 channel analysis (cFos+BLA) with channel C --> BLA
    elif 'NRe' not in  param_list and 'BLA' in param_list:
        df['DAPI'] = data['Num Detections']
        df['cFos'] = data['Num A'] + data['Num AC']
        df['BLA'] = data['Num C'] + data['Num AC']
        df['BLA_cFos'] = data['Num AC']
        df['area'] = data['Area um^2'] / 1e6
    #elif 'Num C' in channels and 'Num B' in channels:        # 3 channel analysis (cFos+BLA+Re) with channel B --> Re and channel C --> BLA
    elif 'NRe' in param_list and 'BLA' in param_list:
        df['DAPI'] = data['Num Detections']
        df['cFos'] = data['Num A'] + data['Num AC'] + data['Num AB']
        df['BLA'] = data['Num C'] + data['Num AC']
        df['BLA_cFos'] = data['Num AC']
        df['NRe'] = data['Num B'] + data['Num AB']
        df['NRe_cFos'] = data['Num AB']
        df['area'] = data['Area um^2'] / 1e6
    
    # Filter the regions that have no detections:
    df = df[df['DAPI']>0]

    return df

#%%
def sum_hemispheres(df, curr_regs):
    '''
    This function sums up the hemispheres of all regions in the dataframe, 
    but only if there are multiple hemispheres for that region.
    
    Input: df = dataframe with all cell counts of a slice, curr_regs = a list with the regions present in the slice.
    Output: hem_df = a dataframe with hemispheres summed.
    '''
    hem_df = pd.DataFrame(np.nan, index=curr_regs, columns=df.columns)
    for region in curr_regs:
        if region+' left' in df.index and region+' right' in df.index:
            sum_hemispheres = df.loc[[region+' left',region+' right']]
        elif region+' left' in df.index and region+' right' not in df.index:
            sum_hemispheres = df.loc[[region+' left']]
        elif region+' left' not in df.index and region+' right' in df.index:
            sum_hemispheres = df.loc[[region+' right']]
        elif region+' left' not in df.index and region+' right' not in df.index:
            sum_hemispheres = df.loc[[region]]
    
        hem_df.loc[region] = sum_hemispheres.sum(axis=0)
    return hem_df

#%%
def import_files(root, combine, param_list):
    """
    This function imports .csv files in the folder specified by 'root', and stores the relevant values into a dictionary called 'mice_dict'.
    Inputs:
        root: path to the folder where all .csv files are located.
        combine: a dictionary with whole regions as keys and a list of subregions as values.
        param_list: list of measured parameters for each region, e.g. ['DAPI','cFos','BLA','area','BLA_cFos'].
    Output:
        mice_dict: a dictionary with all cell counts.
    """
    
    # Make a list of all files in root
    file_list = os.listdir(root)
    # A file '.DS_store' might appear in your list. Remove it:
    if '.DS_Store' in file_list:
        file_list.remove('.DS_Store')
    # Sort the file list
    file_list.sort(key=lambda f: int(''.join(filter(str.isdigit, f))))
    
    # Find the mice that have a file in 'root'
    curr_mouse_list = find_current_mice(file_list)
    
    # Initialize mice_dict
    m_dict = init_dict(curr_mouse_list)
    
    # Loop through csv files in folder:
    for file in file_list:
        print(file)
        mouse_name = file.split('_')[0]
        
        # load data
        data = pd.read_csv(os.path.join(root, file), encoding = 'latin1', index_col = "Name", delimiter=',')
        
        # store raw data in a temporary dataframe
        temp_df = csv_to_dataframe(data, param_list)
        
        # Find the regions and whole regions present in the current slice
        [curr_regs, curr_whole_regs] = find_regs_in_current_slice(temp_df, combine)
        
        # sum hemispheres
        hem_df = sum_hemispheres(temp_df, curr_regs)
    
        # Make an empty table with rows = whole regions in current slice, columns = param_list
        df = pd.DataFrame(np.nan, index=curr_whole_regs, columns=param_list)
        
        # To fill df (=dataframe), subregions are summed
        for curr_whole_reg, curr_subregs in curr_whole_regs.items():
            sum_regs = hem_df.loc[curr_subregs]
            df.loc[curr_whole_reg] = sum_regs.sum(axis=0)
    
        m_dict[mouse_name].append(df)    
        
    return m_dict

#%%
def calculate_ratio(num, den, m_dict, full_mouse_list, whole_regs, thres):
    '''
    This function calculates the average ratio between num (=numerator) and den (=denominator), so ratio = num/den.
    The average is taken over all brain slices in which a brain region appears and which have den above a threshold thres.
    The SEM (variance) is calculated and outputted.
    '''
    # Create empty dataframes for ratios and SEMs.
    ratio_df = pd.DataFrame(np.nan, index=whole_regs, columns=full_mouse_list)
    SEM_df =  pd.DataFrame(np.nan, index=whole_regs, columns=full_mouse_list)
    
    for mouse, df_list in m_dict.items():
        if len(df_list) != 0:
            # Combine slices into a single dataframe
            temp_df = pd.concat([df[[num, den]] for df in df_list])
            # Filter the slices where den is smaller than the threshold
            temp_df = temp_df[temp_df[den] > thres]

            # Fill a dataframe with the ratios
            ratios = pd.DataFrame(np.nan, index=temp_df.index, columns=['ratio'])
            ratios['ratio'] = temp_df[num] / temp_df[den]

            # Take the mean ratio of the slices for each region, and calculate the SEM.
            grouped_ratios = ratios.groupby(ratios.index)
            means = grouped_ratios.mean()
            SEM = grouped_ratios.std(ddof=0) / len(df_list)
            
            # Store ratios and SEMs in the large dataframe
            for region in ratios.index:
                ratio_df.at[region, mouse] = means.loc[region]
                SEM_df.at[region, mouse] = SEM.loc[region]
    
    return ratio_df, SEM_df

#%%
def calculate_chance_ratio(df, tracer):
    '''
    Calculates chance ratio of one slice (dataframe) and stores it in a dataframe 'ratios'.
    '''
    doublepos = tracer + '_cFos'
    ratios = pd.DataFrame(np.nan, index=df.index, columns=['ratio'])
    ratios['ratio'] = (df[doublepos] / df['DAPI']) / ((df['cFos']*df[tracer]) / (df['DAPI']**2))
    return ratios
    
def get_chance_ratio(tracer, m_dict, full_mouse_list, whole_regs, thres):
    '''
    This function calculates the average chance ratio of the 'tracer'. 
    Example: if the tracer is 'BLA', the chance_ratio = ([BLA_cFos]/[DAPI]) / ([cFos]*[BLA])/([DAPI]^2).
    
    The average is taken over all brain slices in which a brain region appears and which have tracer count above a threshold thres.
    The SEM (variance) is calculated and outputted.
    '''
    doublepos = tracer + '_cFos'

    # Create empty dataframes for ratios and SEMs.
    chance_df = pd.DataFrame(np.nan, index=whole_regs, columns=full_mouse_list)
    SEM_df =  pd.DataFrame(np.nan, index=whole_regs, columns=full_mouse_list)

    for mouse, df_list in m_dict.items():
        if len(df_list) != 0:
            # Combine slices into a single dataframe
            temp_df = pd.concat([df[[tracer, doublepos, 'cFos', 'DAPI']] for df in df_list])
            # Filter the slices where the tracer count is smaller than the threshold
            df = temp_df[temp_df[tracer] >= thres]
            
            # Calculate the chance ratios
            chance = calculate_chance_ratio(df, tracer)

            # Take the mean chance of the slices for each region, and calculate the SEM.
            grouped_chance = chance.groupby(chance.index)
            means = grouped_chance.mean()
            SEM = grouped_chance.std(ddof=0) / len(df_list)

            # Store ratios and SEMs in the large dataframe
            for region in means.index:
                chance_df.at[region, mouse] = means.loc[region]
                SEM_df.at[region, mouse] = SEM.loc[region]

    return chance_df, SEM_df
#%%  
def write_results_to_csv(output_directory, output_file, title, recall, extinction, control, tracer,
                         chance_df, chance_SEM_df, 
                         cFos_ratio_df, cFos_SEM_df, 
                         Doublepos_ratio_df, Doublepos_SEM_df):
    '''
    Writes results to an output .csv file.
    '''  
    if not os.path.isdir(output_directory):
        os.mkdir(output_directory)
      
    output_file = os.path.join(output_directory, output_file)
    
    if os.path.isfile(output_file):
        print('WARNING: Output file "{f}" already existed and is overwritten.'.format(f=output_file))
    
    with open(output_file,'w') as f:
        wr = csv.writer(f, quoting=csv.QUOTE_ALL, delimiter=',')
        wr.writerow([title])
        wr.writerow('\n')
        wr.writerow(['Recall:']+recall)
        wr.writerow(['Extinction:']+extinction)
        wr.writerow(['Control:']+control)

    with open(output_file,'a+') as f:
        f.write('\n {i} chance ratio \n'.format(i=tracer))
    chance_df.to_csv(output_file, mode='a')
    
    with open(output_file,'a+') as f:
        f.write('\n {i} chance ratio SEM \n'.format(i=tracer))
    chance_SEM_df.to_csv(output_file, mode='a')
        
    with open(output_file,'a+') as f:
        f.write('\n cFos ratio \n')
    cFos_ratio_df.to_csv(output_file, mode='a')
    
    with open(output_file,'a+') as f:
        f.write('\n cFos ratio SEM \n')
    cFos_SEM_df.to_csv(output_file, mode='a')

    with open(output_file,'a+') as f:
        f.write('\n {i}_cFos/{i} (doublepos) ratio \n'.format(i=tracer))
    Doublepos_ratio_df.to_csv(output_file, mode='a')
    
    with open(output_file,'a+') as f:
        f.write('\n {i}_cFos/{i} (doublepos) SEM \n'.format(i=tracer))
    Doublepos_SEM_df.to_csv(output_file, mode='a')
    
    f.close()
#%%
def output_raw_numbers(output_directory, m_dict, tracer, region_name):
    '''
    Outputs raw numbers (the numbers used for caclulating the chance ratio).
    '''
    
    if not os.path.isdir(output_directory):
        os.mkdir(output_directory)
    
    raw_number_directory = os.path.join(output_directory, 'Raw_numbers'+'_'+region_name+'_'+tracer)
    
    if os.path.isdir(raw_number_directory):
        shutil.rmtree(raw_number_directory)
        print('WARNING: directory "{d}" already existed and is overwritten'.format(d=raw_number_directory))
        
    os.mkdir(raw_number_directory)
    
    for mouse, df_list in m_dict.items():
        file_name = os.path.join(raw_number_directory, mouse + '.csv')
        
        with open(file_name,'w') as f:
            wr = csv.writer(f, quoting=csv.QUOTE_ALL, delimiter=',')
            wr.writerow([mouse])
            wr.writerow('\n')
        
        s = 0
        for df in df_list:
            s = s + 1
            df = df.drop(columns='area')
            with open(file_name,'a+') as f:
                f.write('\n slice {i} \n'.format(i=s))
            df.to_csv(file_name, mode='a', header=True)
    
    
    