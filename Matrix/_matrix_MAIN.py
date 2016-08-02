#!/usr/bin/python
# -*- coding: utf-8 -*-


"""
Created on Mon Jan 18 09:44:56 2016

@author: stefano
"""


#++++++++++++++ Requested Package(s) Import ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
import os

import matrix_configure_module

from matrix_preprocessing_assoFile_module import loadAssoFile

from matrix_preprocessing_dataSources_module import getLaunchPathDict
from matrix_preprocessing_dataLoading_module import loadData

from matrix_processing_computeMatrixes_module import buildSeqCountMatrix, buildShsCountMatrix, buildBarcodeCountMatrix, buildCellCountMatrix, buildFragmentEstimateMatrix

from matrix_output_module import buildOutputPath, writeMatrix
#++++++++++++++++++++++ Global Vars from globModule ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

### Import Functions
humanSorted = matrix_configure_module.humanSorted

### Screen print configs 
verbose = matrix_configure_module.verbose
verbosePrint = matrix_configure_module.verbosePrint

### Association File Loading configs - assoModule
asso_folder = matrix_configure_module.asso_folder
asso_file_name = matrix_configure_module.asso_file_name
asso_delimiter = matrix_configure_module.asso_delimiter

### Data Loading configs - dataModule
dataset_tuple_list = matrix_configure_module.dataset_tuple_list
drop_headers = matrix_configure_module.drop_headers
compression = matrix_configure_module.compression

### COMMON OUTPUT GROUND DIR
common_output_ground_dir = matrix_configure_module.common_output_ground_dir

### Matrix output configs - outputModule
matrix_outfolder = matrix_configure_module.matrix_outfolder
matrix_files_delimiter = matrix_configure_module.matrix_files_delimiter
relabelling = matrix_configure_module.relabelling  # Boolean


#++++++++++++++++++++++++ CODE +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#


### Load Association File Data #######################################
asso_dict = loadAssoFile(asso_file_name, asso_folder, asso_delimiter)
######################################################################

### Load Data ########################################################
launch_path_dict = getLaunchPathDict(dataset_tuple_list)
any_df = loadData(launch_path_dict, drop_headers, compression)
######################################################################

### Compute matrixes ############################################################################################################
seqCount_matrix, ShsCount_matrix, barcodeCount_matrix, cellCount_matrix, fragmentEstimate_matrix = None, None, None, None, None
verbosePrint("\n>>> Computing matrixes:")
verbosePrint("> seqCount matrix ...")
seqCount_matrix = buildSeqCountMatrix(any_df)
verbosePrint("> ShsCount matrix ...")
ShsCount_matrix = buildShsCountMatrix(any_df)
verbosePrint("> barcodeCount matrix ...")
barcodeCount_matrix = buildBarcodeCountMatrix(any_df)
verbosePrint("> cellCount matrix ...")
cellCount_matrix = buildCellCountMatrix(any_df)
verbosePrint("> fragmentEstimate matrix ...")
fragmentEstimate_matrix = buildFragmentEstimateMatrix(any_df)
verbosePrint(">>> Done!")
#################################################################################################################################

### Output ################################################################################################################################################################
metadata = None  # No relabelling
if relabelling:
    metadata = asso_dict  # relabelling columns during output generation

verbosePrint("\n>>> Export matrixes ...")
# create a sub folder of common_output_ground_dir where write matrixes
matrix_files_outdir = buildOutputPath(common_output_ground_dir, matrix_outfolder)
seqCount_matrix_outPath, ShsCount_matrix_outPath, barcodeCount_matrix_outPath, cellCount_matrix_outPath, fragmentEstimate_matrix_outPath  = None, None, None, None, None
# seqCount matrix
seqCount_matrix_outPath = writeMatrix(seqCount_matrix, os.path.join(matrix_files_outdir, "seqCount_matrix.tsv"), matrix_files_delimiter, metadata=metadata)
# ShsCount matrix
ShsCount_matrix_outPath = writeMatrix(ShsCount_matrix, os.path.join(matrix_files_outdir, "ShsCount_matrix.tsv"), matrix_files_delimiter, metadata=metadata)
# barcodeCount matrix
barcodeCount_matrix_outPath = writeMatrix(barcodeCount_matrix, os.path.join(matrix_files_outdir, "barcodeCount_matrix.tsv"), matrix_files_delimiter, metadata=metadata)
# cellCount matrix
cellCount_matrix_outPath = writeMatrix(cellCount_matrix, os.path.join(matrix_files_outdir, "cellCount_matrix.tsv"), matrix_files_delimiter, metadata=metadata)
# fragmentEstimate matrix
fragmentEstimate_matrix_outPath = writeMatrix(fragmentEstimate_matrix, os.path.join(matrix_files_outdir, "fragmentEstimate_matrix.tsv"), matrix_files_delimiter, metadata=metadata)
verbosePrint(">>> Matrix Files Created!")

# Note for normalization
#seqCount_matrix_norm = seqCount_matrix.apply(lambda x: x/x.sum())
#ShsCount_matrix_norm = ShsCount_matrix.apply(lambda x: x/x.sum())
#barcodeCount_matrix_norm = barcodeCount_matrix.apply(lambda x: x/x.sum())
#cellCount_matrix_norm = cellCount_matrix.apply(lambda x: x/x.sum())
#fragmentEstimate_matrix_norm = fragmentEstimate_matrix.apply(lambda x: x/x.sum())

##########################################################################################################################################################################


#+++++++++++++++++++++++++ END CODE ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#






#++++++++++++++++++++++++ LOCAL TEST +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#



#++++++++++++++++++++++++ DEVEL ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#


