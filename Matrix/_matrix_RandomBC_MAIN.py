# -*- coding: utf-8 -*-


"""
Created on Mon Jan 18 09:44:56 2016

@author: stefano
"""


#++++++++++++++ Requested Package(s) Import +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

import os

import matrix_RandomBC_globModule
import matrix_RandomBC_assoModule
import matrix_RandomBC_dataModule


#++++++++++++++++++++++ Global Vars ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

# Screen print
verbose = matrix_RandomBC_globModule.verbose
verbosePrint = matrix_RandomBC_globModule.verbosePrint

# Association File - assoModule
asso_folder = matrix_RandomBC_globModule.asso_folder
asso_file_name = matrix_RandomBC_globModule.asso_file_name
asso_delimiter = matrix_RandomBC_globModule.asso_delimiter

# Data - dataModule
ground_dir = matrix_RandomBC_globModule.ground_dir
DISEASE = matrix_RandomBC_globModule.DISEASE
PATIENT = matrix_RandomBC_globModule.PATIENT
POOL = matrix_RandomBC_globModule.POOL
data_files_delimiter = matrix_RandomBC_globModule.data_files_delimiter
data_files_name_filter = matrix_RandomBC_globModule.data_files_name_filter


#++++++++++++++++++++++++ CODE +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

# OVERRIDE ASSO VARS FOR LOCAL TEST
asso_folder = "/home/stefano/Desktop/RandomBC_matrix_development/test_input/asso"
asso_file_name = "asso.assayvalidation.lane1.tsv"
asso_delimiter = '\t'
### Load Association File Data #################################################################
asso_dict = matrix_RandomBC_assoModule.loadAssoFile(asso_file_name, asso_folder, asso_delimiter)
################################################################################################

# OVERRIDE DATA VARS FOR LOCAL TEST
data_files_delimiter = '\t'
data_files_name_filter = ".randomBC.tsv"
develop_input_data_path = "/home/stefano/Desktop/RandomBC_matrix_development/test_input/data"
# TMP CODE FOR LOCAL TEST
filtered_dir_content = matrix_RandomBC_dataModule.listDir(develop_input_data_path, name_filter=data_files_name_filter)
verbosePrint("\n\n>>> Loading data ...")
verbosePrint("> develop_input_data_path: {develop_input_data_path}".format(develop_input_data_path=str(develop_input_data_path)))
verbosePrint("> exploited substring for data detection: '{data_files_name_filter}'".format(data_files_name_filter=str(data_files_name_filter)))
verbosePrint("> n data files detected: {n_files}".format(n_files=str(len(filtered_dir_content))))
verbosePrint("> data file list: {filtered_dir_content}".format(filtered_dir_content=str(filtered_dir_content)))
verbosePrint("")
POOL_IS_dict = {}
POOL_alldata_dict = {}
for path in filtered_dir_content:
    filename = str(os.path.basename(path))
    barcode = ".".join((filename.split("."))[:2])
    verbosePrint("> Processing {filename}, barcode={barcode} ...".format(filename=str(filename), barcode=str(barcode)))
    data_file_nested_list = matrix_RandomBC_dataModule.loadFile (path, data_files_delimiter)
    alldata_dict, IS_dict = matrix_RandomBC_dataModule.arrangeData(data_file_nested_list)
    POOL_IS_dict[barcode] = IS_dict
    POOL_alldata_dict[barcode] = alldata_dict
verbosePrint("\n>>> Data Loaded!\n")
### Load Data ################################################################################################################################################
# POOL_alldata_dict, POOL_IS_dict = matrix_RandomBC_dataModule.loadDataFiles(ground_dir, DISEASE, PATIENT, POOL, data_files_name_filter, data_files_delimiter)
##############################################################################################################################################################




#++++++++++++++++++++++++ DEVEL +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

# import packages
import pandas as pd
import numpy as np
import collections

def flatten(d, parent_key='', sep='@'):
    items = []
    for k, v in d.items():
        new_key = parent_key + sep + str(k) if parent_key else k
        if isinstance(v, collections.MutableMapping):
            items.extend(flatten(v, new_key, sep=sep).items())
        else:
            items.append((tuple(new_key.split(sep)), v))
    return dict(items)

# flatten data
flat_POOL_IS_dict = flatten(POOL_IS_dict)
l=[]
[l.append(list(k) + [v]) for k,v in flat_POOL_IS_dict.items()]

# load as dataframe
df_flat = pd.DataFrame(l, columns=['barcode', 'genomic_coordinates', 'shearsite', 'randomBC', 'seq_count'])

# pivot to get sequence counts of ISs
df_flat_seqCount = pd.pivot_table(df_flat, index='genomic_coordinates', columns='barcode', values='seq_count', aggfunc='sum', fill_value=0)

# pivot to get counts of distinct SHEARSITE-RANDOMTAG couples for each ISs
df_flat_cellCount = pd.pivot_table(df_flat, index='genomic_coordinates', columns='barcode', values='seq_count', aggfunc=np.count_nonzero, fill_value=0)

