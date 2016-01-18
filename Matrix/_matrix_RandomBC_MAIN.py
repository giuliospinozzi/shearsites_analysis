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
import matrix_RandomBC_processingModule
import matrix_RandomBC_outputModule


#++++++++++++++++++++++ Global Vars from globModule +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

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

# Output - outputModule
#ground_dir, DISEASE, PATIENT, POOL as for Data
outfolder = matrix_RandomBC_globModule.outfolder
out_files_delimiter = matrix_RandomBC_globModule.out_files_delimiter


#++++++++++++++++++++++++ CODE +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#


#### Load Association File Data ###################################################################
#asso_dict = matrix_RandomBC_assoModule.loadAssoFile(asso_file_name, asso_folder, asso_delimiter)
###################################################################################################
#
#### Load Data ################################################################################################################################################
#POOL_alldata_dict, POOL_IS_dict = matrix_RandomBC_dataModule.loadDataFiles(ground_dir, DISEASE, PATIENT, POOL, data_files_name_filter, data_files_delimiter)
###############################################################################################################################################################
#
#### Process Data ###########################################################
#df = matrix_RandomBC_processingModule.buildDataFrame(POOL_IS_dict)
#seqCount_matrix = matrix_RandomBC_processingModule.buildSeqCountMatrix(df)
#cellCount_matrix = matrix_RandomBC_processingModule.buildCellCountMatrix(df)
#ShsCount_matrix = matrix_RandomBC_processingModule.buildShsCountMatrix(df)
#############################################################################
#
#### Output ##############################################################################################
#OUTDIR = matrix_RandomBC_outputModule.buildOutputPath(ground_dir, DISEASE, PATIENT, POOL, outfolder)
#seqCount_matrix_outPath = os.path.normpath(os.path.join(OUTDIR, "seqCount_matrix.tsv"))
#matrix_RandomBC_outputModule.writeMatrix(seqCount_matrix, seqCount_matrix_outPath, out_files_delimiter)
#cellCount_matrix_outPath = os.path.normpath(os.path.join(OUTDIR, "cellCount_matrix.tsv"))
#matrix_RandomBC_outputModule.writeMatrix(cellCount_matrix, cellCount_matrix_outPath, out_files_delimiter)
#ShsCount_matrix_outPath = os.path.normpath(os.path.join(OUTDIR, "ShsCount_matrix.tsv"))
#matrix_RandomBC_outputModule.writeMatrix(ShsCount_matrix, ShsCount_matrix_outPath, out_files_delimiter)
##########################################################################################################

#++++++++++++++++++++++++ LOCAL TEST +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

### Load Association File Data
# OVERRIDE ASSO VARS FOR LOCAL TEST
asso_folder = "/home/stefano/Desktop/RandomBC_matrix_development/test_input/asso"
asso_file_name = "asso.assayvalidation.lane1.tsv"
asso_delimiter = '\t'
# load
asso_dict = matrix_RandomBC_assoModule.loadAssoFile(asso_file_name, asso_folder, asso_delimiter)


### Load Data
# OVERRIDE DATA VARS FOR LOCAL TEST
data_files_delimiter = '\t'
data_files_name_filter = ".randomBC.tsv"
develop_input_data_path = "/home/stefano/Desktop/RandomBC_matrix_development/test_input/data"
# TMP CODE FOR LOCAL TEST
# in place of POOL_alldata_dict, POOL_IS_dict = matrix_RandomBC_dataModule.loadDataFiles(ground_dir, DISEASE, PATIENT, POOL, data_files_name_filter, data_files_delimiter)
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

#++++++++++++++++++++++++ DEVEL ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

### Process Data
df = matrix_RandomBC_processingModule.buildDataFrame(POOL_IS_dict)
# df[['shearsite']] = df[['shearsite']].astype(int) # if cast is needed
seqCount_matrix = matrix_RandomBC_processingModule.buildSeqCountMatrix(df)
cellCount_matrix = matrix_RandomBC_processingModule.buildCellCountMatrix(df)
ShsCount_matrix = matrix_RandomBC_processingModule.buildShsCountMatrix(df)

### Output
#OUTDIR = matrix_RandomBC_outputModule.buildOutputPath(ground_dir, DISEASE, PATIENT, POOL, outfolder)
OUTDIR = "/home/stefano/Desktop/RandomBC_matrix_development/test_output"

seqCount_matrix_outPath = os.path.normpath(os.path.join(OUTDIR, "seqCount_matrix.tsv"))
matrix_RandomBC_outputModule.writeMatrix(seqCount_matrix, seqCount_matrix_outPath, out_files_delimiter)

cellCount_matrix_outPath = os.path.normpath(os.path.join(OUTDIR, "cellCount_matrix.tsv"))
matrix_RandomBC_outputModule.writeMatrix(cellCount_matrix, cellCount_matrix_outPath, out_files_delimiter)

ShsCount_matrix_outPath = os.path.normpath(os.path.join(OUTDIR, "ShsCount_matrix.tsv"))
matrix_RandomBC_outputModule.writeMatrix(ShsCount_matrix, ShsCount_matrix_outPath, out_files_delimiter)









