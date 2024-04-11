#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 20:33:46 2023

@author: zmk214
"""

# add options to select file withWord, withoutWord in readInProfiles

# libraries
import os
import pandas as pd
import numpy as np
from sklearn.utils import shuffle

def readInProfiles(path, setName="", withAug=False, withWord=False, withoutWord=False):
    
    # input: path of the folder (note that all files have to be .csv files with header)
    # process: all files in the folder will be read and concat together
    # output: pandas dataframe object contains all data from files in the folder
    
    print(setName)
    
    if withWord == False:
        dirlsnames = os.listdir(path)
    elif withWord == True:
        print("not correct param")
    else:
        dirlsnames = os.listdir(path)
        dirlsnames = [name for name in dirlsnames if all(substring in name for substring in withWord)]

        
    if withoutWord == False:
        print("no remove word")
    elif withoutWord == True:
        print("not correct param")
    else:
        dirlsnames = [name for name in dirlsnames if all(substring not in name for substring in withoutWord)]
    
    if '.DS_Store' in dirlsnames:
        dirlsnames.remove('.DS_Store')

    names_with_aug = [name for name in dirlsnames if 'aug' in name]
    names_without_aug = [name for name in dirlsnames if 'aug' not in name]
    
    if withAug == True:
        dirls = names_with_aug + names_without_aug
    else:
        dirls = names_without_aug
         
    dirls.sort()
    
    data_in = []
    data_in = pd.DataFrame(data_in)
    
    for i in dirls:
        read_in = pd.read_csv(path + i)
        print(i, "contains", len(read_in), "rows")
        data_in = pd.concat([data_in, read_in])
    
    print("In total, the dataset contains", len(data_in), "rows", "\n")
    return data_in

def readInMetadata(path, setName="", withAug=False):
    
    # all the same as readInProfiles() for now
    # input: path of the folder (note that all files have to be .csv files with header)
    # process: all files in the folder will be read and concat together
    # output: pandas dataframe object contains all data from files in the folder
    
    print(setName)
    
    dirlsnames = os.listdir(path)
    
    if '.DS_Store' in dirlsnames:
        dirlsnames.remove('.DS_Store')
    
    names_with_aug = [name for name in dirlsnames if 'aug' in name]
    names_without_aug = [name for name in dirlsnames if 'aug' not in name]
    
    if withAug == True:
        dirls = names_with_aug + names_without_aug
    else:
        dirls = names_without_aug
         
    dirls.sort()
    
    data_in = []
    data_in = pd.DataFrame(data_in)
    
    for i in dirls:
        read_in = pd.read_csv(path + i)
        print(i, "contains", len(read_in), "rows")
        data_in = pd.concat([data_in, read_in])
    
    print("In total, the dataset contains", len(data_in), "rows", "\n")
    return data_in

def generateModelInput(pos, neg, equalNumber=False, addShuffle=True):
    
    # input: dataframe (output from readInProfiles)
    # process: 1) combine pos and neg sets 
    #          2) if equalNumber=True, randomly select the larger set to have equal number with another set
    #          3) create the label (y_) 
    # output: numpy array objects of x=data and y=label already shuffle (default)
    
    if equalNumber==True:
        if len(pos) > len(neg):
            pos = pos.sample(n=len(neg))
        elif len(neg) > len(pos):
            neg = neg.sample(n=len(pos))
            
    x_df = pd.concat([pos, neg])
    x_np = x_df.to_numpy() 
    y_np = np.array([1]*len(pos) + [0]*len(neg))
    
    if addShuffle==True:
        x_np, y_np = shuffle(x_np, y_np)
        
    return x_np, y_np

def generateModelInputWithMetadata(pos, neg, posmt, negmt, equalNumber=False, addShuffle=True):
    
    # input: dataframe (output from readInProfiles)
    # process: 1) combine pos and neg sets 
    #          2) if equalNumber=True, randomly select the larger set to have equal number with another set
    #          3) create the label (y_) 
    # output: numpy array objects of x=data and y=label already shuffle (default)
    
    if equalNumber==True:
        if len(pos) > len(neg):
            pos = pos.sample(n=len(neg))
            posmt = posmt.loc[pos.index]
            pos.sort_index(inplace=True)
            posmt.sort_index(inplace=True)
        elif len(neg) > len(pos):
            neg = neg.sample(n=len(pos))
            negmt = negmt.loc[neg.index]
            neg.sort_index(inplace=True)
            negmt.sort_index(inplace=True)
    
    x_df = pd.concat([pos, neg])
    x_np = x_df.to_numpy() 
    y_np = np.array([1]*len(pos) + [0]*len(neg))
    
    input_df_dict = {"positive_prfl": pos,
                     "positive_mtdt": posmt, 
                     "negative_prfl": neg,
                     "negative_mtdt": negmt} 
    
    if addShuffle==True:
        x_np, y_np = shuffle(x_np, y_np)
        
    return x_np, y_np, input_df_dict

def seperateAnnotation(prfl, mtdt, txtype, val_frac=0.2):
    
    new_mtdt_index = mtdt.index[mtdt['atac_txType'] == txtype].tolist()
    new_mtdt = mtdt.loc[new_mtdt_index, ]
    new_prfl = prfl.loc[new_mtdt_index, ]
    print(txtype, len(new_mtdt_index))
    
    new_prfl_val = new_prfl.sample(frac=val_frac, replace=False, random_state=214)
    new_mtdt_val = new_mtdt[ new_mtdt.index.isin(new_prfl_val.index) ]
    
    new_prfl_tra = new_prfl[ ~new_prfl.index.isin(new_prfl_val.index) ]
    new_mtdt_tra = new_mtdt[ ~new_mtdt.index.isin(new_prfl_val.index) ]
    print("new_prfl_tra", len(new_prfl_tra))
    print("new_prfl_val", len(new_prfl_val))
    
    return new_prfl, new_prfl_tra, new_prfl_val, new_mtdt, new_mtdt_tra, new_mtdt_val