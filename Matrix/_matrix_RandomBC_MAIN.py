#!/usr/bin/python
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
import matrix_RandomBC_filterModule
import matrix_RandomBC_outputModule


#++++++++++++++++++++++ Global Vars from globModule +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

### Import Functions
humanSorted = matrix_RandomBC_globModule.humanSorted

### Screen print configs 
verbose = matrix_RandomBC_globModule.verbose
verbosePrint = matrix_RandomBC_globModule.verbosePrint

### Association File Loading configs - assoModule
asso_folder = matrix_RandomBC_globModule.asso_folder
asso_file_name = matrix_RandomBC_globModule.asso_file_name
asso_delimiter = matrix_RandomBC_globModule.asso_delimiter

### Data Loading configs - dataModule
ground_dir = matrix_RandomBC_globModule.ground_dir
DISEASE = matrix_RandomBC_globModule.DISEASE
PATIENT = matrix_RandomBC_globModule.PATIENT
POOL = matrix_RandomBC_globModule.POOL
data_files_delimiter = matrix_RandomBC_globModule.data_files_delimiter
data_files_name_filter = matrix_RandomBC_globModule.data_files_name_filter

### Cleaning Dataframe configs - filterModule
clean_df = matrix_RandomBC_globModule.clean_df
ED_rule = matrix_RandomBC_globModule.ED_rule
ED_inside_ShS = matrix_RandomBC_globModule.ED_inside_ShS

### COMMON OUTPUT GROUND DIR
common_output_ground_dir = matrix_RandomBC_globModule.common_output_ground_dir

### Matrix output configs - outputModule
matrix_outfolder = matrix_RandomBC_globModule.matrix_outfolder
matrix_files_delimiter = matrix_RandomBC_globModule.matrix_files_delimiter
relabelling = matrix_RandomBC_globModule.relabelling  # Boolean


#++++++++++++++++++++++++ CODE +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#


### Load Association File Data ###################################################################
asso_dict = matrix_RandomBC_assoModule.loadAssoFile(asso_file_name, asso_folder, asso_delimiter)
##################################################################################################

### Load Data ################################################################################################################################################
POOL_alldata_dict, POOL_IS_dict = matrix_RandomBC_dataModule.loadDataFiles(ground_dir, DISEASE, PATIENT, POOL, data_files_name_filter, data_files_delimiter)
##############################################################################################################################################################

### Process Data ##################################################################
verbosePrint("\n>>> Shape data as DataFrame ...")
df = matrix_RandomBC_processingModule.buildDataFrame(POOL_IS_dict)
#df = matrix_RandomBC_processingModule.buildExhaustiveDataFrame(POOL_alldata_dict)
verbosePrint(">>> Dataframe built!")
###################################################################################

### Filter Data ####################################################################################################
if clean_df:
    verbosePrint("\n>>> Cleaning DataFrame ...")
    df = matrix_RandomBC_filterModule.filterBy_randomBC_EditDistance(df, inside_ShS=ED_inside_ShS, ED_rule=ED_rule)
    # Here new filters!
    verbosePrint(">>> Done!")
####################################################################################################################

### Compute matrixes ############################################################################################################
seqCount_matrix, ShsCount_matrix, barcodeCount_matrix, cellCount_matrix, fragmentEstimate_matrix = None, None, None, None, None
verbosePrint("\n>>> Computing matrixes:")
verbosePrint("> seqCount matrix ...")
seqCount_matrix = matrix_RandomBC_processingModule.buildSeqCountMatrix(df)
verbosePrint("> ShsCount matrix ...")
ShsCount_matrix = matrix_RandomBC_processingModule.buildShsCountMatrix(df)
verbosePrint("> barcodeCount matrix ...")
barcodeCount_matrix = matrix_RandomBC_processingModule.buildBarcodeCountMatrix(df)
verbosePrint("> cellCount matrix ...")
cellCount_matrix = matrix_RandomBC_processingModule.buildCellCountMatrix(df)
verbosePrint("> fragmentEstimate matrix ...")
fragmentEstimate_matrix = matrix_RandomBC_processingModule.buildFragmentEstimateMatrix(df)
verbosePrint(">>> Done!")
#################################################################################################################################

### Output ###################################################################################################################################################
metadata = None  # No relabelling
if relabelling:
    metadata = asso_dict  # relabelling columns during output generation
verbosePrint("\n>>> Export matrixes ...")
# Export: make OUTDIR -> make outPath -> call writeMatrix
OUTDIR = matrix_RandomBC_outputModule.buildOutputPath(common_output_ground_dir, matrix_outfolder)
seqCount_matrix_outPath = os.path.normpath(os.path.join(OUTDIR, "seqCount_matrix.tsv"))
matrix_RandomBC_outputModule.writeMatrix(seqCount_matrix, seqCount_matrix_outPath, matrix_files_delimiter, metadata=metadata)
ShsCount_matrix_outPath = os.path.normpath(os.path.join(OUTDIR, "ShsCount_matrix.tsv"))
matrix_RandomBC_outputModule.writeMatrix(ShsCount_matrix, ShsCount_matrix_outPath, matrix_files_delimiter, metadata=metadata)
barcodeCount_matrix_outPath = os.path.normpath(os.path.join(OUTDIR, "barcodeCount_matrix.tsv"))
matrix_RandomBC_outputModule.writeMatrix(barcodeCount_matrix, barcodeCount_matrix_outPath, matrix_files_delimiter, metadata=metadata)
cellCount_matrix_outPath = os.path.normpath(os.path.join(OUTDIR, "cellCount_matrix.tsv"))
matrix_RandomBC_outputModule.writeMatrix(cellCount_matrix, cellCount_matrix_outPath, matrix_files_delimiter, metadata=metadata)
fragmentEstimate_matrix_outPath = os.path.normpath(os.path.join(OUTDIR, "fragmentEstimate_matrix.tsv"))
matrix_RandomBC_outputModule.writeMatrix(fragmentEstimate_matrix, fragmentEstimate_matrix_outPath, matrix_files_delimiter, metadata=metadata)
verbosePrint(">>> Matrix Files Created!")

# Note for normalization
#seqCount_matrix_norm = seqCount_matrix.apply(lambda x: x/x.sum())
#ShsCount_matrix_norm = ShsCount_matrix.apply(lambda x: x/x.sum())
#barcodeCount_matrix_norm = barcodeCount_matrix.apply(lambda x: x/x.sum())
#cellCount_matrix_norm = cellCount_matrix.apply(lambda x: x/x.sum())
#fragmentEstimate_matrix_norm = fragmentEstimate_matrix.apply(lambda x: x/x.sum())

#############################################################################################################################################################


#+++++++++++++++++++++++++ END CODE ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#






##++++++++++++++++++++++++ LOCAL TEST +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#
#
#### Load Association File Data ##################################################################
## OVERRIDE ASSO VARS FOR LOCAL TEST
#asso_folder = "/home/stefano/Desktop/RandomBC_matrix_development/test_input/asso"
#asso_file_name = "asso.assayvalidation.lane1.tsv"
#asso_delimiter = '\t'
## load
#asso_dict = matrix_RandomBC_assoModule.loadAssoFile(asso_file_name, asso_folder, asso_delimiter)
##################################################################################################
#
#### Load Data ###############################################################################################################################################################
## OVERRIDE DATA VARS FOR LOCAL TEST
#data_files_delimiter = '\t'
#data_files_name_filter = ".randomBC.tsv"
#develop_input_data_path = "/home/stefano/Desktop/RandomBC_matrix_development/test_input/data"
#
## in place of POOL_alldata_dict, POOL_IS_dict = matrix_RandomBC_dataModule.loadDataFiles(ground_dir, DISEASE, PATIENT, POOL, data_files_name_filter, data_files_delimiter)
#filtered_dir_content = matrix_RandomBC_dataModule.listDir(develop_input_data_path, name_filter=data_files_name_filter)
#verbosePrint("\n\n>>> Loading data ...")
#verbosePrint("> develop_input_data_path: {develop_input_data_path}".format(develop_input_data_path=str(develop_input_data_path)))
#verbosePrint("> exploited substring for data detection: '{data_files_name_filter}'".format(data_files_name_filter=str(data_files_name_filter)))
#verbosePrint("> n data files detected: {n_files}".format(n_files=str(len(filtered_dir_content))))
#verbosePrint("> data file list: {filtered_dir_content}".format(filtered_dir_content=str(filtered_dir_content)))
#verbosePrint("")
#POOL_IS_dict = {}
#POOL_alldata_dict = {}
#for path in filtered_dir_content:
#    filename = str(os.path.basename(path))
#    barcode = ".".join((filename.split("."))[:2])
#    verbosePrint("> Processing {filename}, barcode={barcode} ...".format(filename=str(filename), barcode=str(barcode)))
#    data_file_nested_list = matrix_RandomBC_dataModule.loadFile (path, data_files_delimiter)
#    alldata_dict, IS_dict = matrix_RandomBC_dataModule.arrangeData(data_file_nested_list)
#    POOL_IS_dict[barcode] = IS_dict
#    POOL_alldata_dict[barcode] = alldata_dict
#verbosePrint("\n>>> Data Loaded!\n")
##############################################################################################################################################################################
#
#### Process Data ##################################################################
#verbosePrint("\n>>> Shape data as DataFrame ...")
#df = matrix_RandomBC_processingModule.buildDataFrame(POOL_IS_dict)
##df = matrix_RandomBC_processingModule.buildExhaustiveDataFrame(POOL_alldata_dict)
#verbosePrint(">>> Dataframe built!")
####################################################################################
#
#
#
##++++++++++++++++++++++++ DEVEL ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
