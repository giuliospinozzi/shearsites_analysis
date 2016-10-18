#!/usr/bin/python
# -*- coding: utf-8 -*-


"""
Created on Thu Oct 13 14:29:40 2016

@author: Stefano Brasca
"""

__author__ = "Stefano Brasca"
__copyright__ = "SR-TIGET"
__credits__ = ["Stefano Brasca", "Andrea Calabria", "Giulio Spinozzi"]
__version__ = "1.0"
__maintainer__ = "Stefano Brasca"
__email__ = "brasca.stefano@hsr.it"
__status__ = "Testing"


#++++++++++++++ Requested Package(s) Import ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

### Check requirements, import configuration and PARSE ARGS
import matrix_configure_module

### Import Basic Functions
humanSorted = matrix_configure_module.humanSorted
verbosePrint = matrix_configure_module.verbosePrint
from matrix_preprocessing_dataSources_module import getLaunchPathDict
from matrix_preprocessing_dataLoading_module import loadData
from matrix_processing_computeMatrixes_module import buildSeqCountMatrix, buildShsCountMatrix, buildBarcodeCountMatrix, buildCellCountMatrix, buildFragmentEstimateMatrix
from matrix_output_module import buildOutputPath, writeMatrix


#++++++++++++++++++++++ Global Vars from configure_module ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

### Screen print configs 
verbose = matrix_configure_module.verbose

### Data Loading configs - dataModule
dataset_tuple_list = matrix_configure_module.dataset_tuple_list
drop_headers = matrix_configure_module.drop_headers
compression = matrix_configure_module.compression

### ISs computation - matrix_processing_ISsMethods_module
do_ISs = matrix_configure_module.do_ISs
# ensembles config
ensembles_per_sample = matrix_configure_module.ensembles_per_sample
ensembles_max_dist = matrix_configure_module.ensembles_max_dist
ensembles_max_span = matrix_configure_module.ensembles_max_span
# ISs method
ISs_method = matrix_configure_module.ISs_method

### Filter Data configs - matrix_preprocessing_filterData_module
filter_data = matrix_configure_module.filter_data
filter_by_ED = matrix_configure_module.filter_by_ED
inside_ShS = matrix_configure_module.inside_ShS
ED_treshold = matrix_configure_module.ED_treshold

### COMMON OUTPUT GROUND DIR
common_output_ground_dir = matrix_configure_module.common_output_ground_dir

### Matrix output configs - outputModule
matrix_outfolder = matrix_configure_module.matrix_outfolder
matrix_files_delimiter = matrix_configure_module.matrix_files_delimiter


#++++++++++++++++++++++++ CODE +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

verbosePrint("\n[START]")

### Print Config ##########################################################################
verbosePrint("\n>>> Configuring ...")
verbosePrint("    DATA:")
verbosePrint("    * dataset_tuple_list: {x}".format(x=str(dataset_tuple_list)))
verbosePrint("    * drop_headers: {x}".format(x=str(drop_headers)))
verbosePrint("    * compression: {x}".format(x=str(compression)))
verbosePrint("    ISs:")
verbosePrint("    * do_ISs: {x}".format(x=str(do_ISs)))
if do_ISs:
    verbosePrint("     * ensembles_per_sample: {x}".format(x=str(ensembles_per_sample)))
    verbosePrint("     * ensembles_max_dist: {x}".format(x=str(ensembles_max_dist)))
    verbosePrint("     * ensembles_max_span: {x}".format(x=str(ensembles_max_span)))
    verbosePrint("     * ISs_method: {x}".format(x=str(ISs_method)))
verbosePrint("    FILTER:")
verbosePrint("    * filter_data: {x}".format(x=str(filter_data)))
if filter_data:
    verbosePrint("     * filter_by_ED: {x}".format(x=str(filter_by_ED)))
    verbosePrint("     * inside_ShS: {x}".format(x=str(inside_ShS)))
    verbosePrint("     * ED_treshold: {x}".format(x=str(ED_treshold)))
verbosePrint("    OUTPUT:")
verbosePrint("    * common_output_ground_dir: {x}".format(x=str(common_output_ground_dir)))
if matrix_outfolder != '':
    verbosePrint("    * matrix_outfolder: {x}".format(x=str(matrix_outfolder)))
if matrix_files_delimiter == '\t':
    verbosePrint(r'''    * matrix_files_delimiter: \t''')
else:
    verbosePrint("    * matrix_files_delimiter: {x}".format(x=str(matrix_files_delimiter)))
verbosePrint(">>> Done!")
###########################################################################################

### Check (or create) matrix_files_outdir ########################################
verbosePrint("\n>>> Setting up OUTDIR ...")
matrix_files_outdir = buildOutputPath(common_output_ground_dir, matrix_outfolder)
verbosePrint(">>> Done!")
##################################################################################

### Load Data ################################################
launch_path_dict = getLaunchPathDict(dataset_tuple_list)
any_df = loadData(launch_path_dict, drop_headers, compression)
##############################################################

### Compute ISs ################################################################################################################################################################
if do_ISs:
    verbosePrint("\n>>> Computing ISs ...")
    from matrix_processing_ISsMethods_module import compute_ISs
    any_df = compute_ISs(any_df, ensembles_per_sample=ensembles_per_sample, ensembles_max_dist=ensembles_max_dist, ensembles_max_span=ensembles_max_span, ISs_method=ISs_method)
    verbosePrint(">>> Done!")
################################################################################################################################################################################

### Filter Data #######################################################################################################################################
if filter_data:
    verbosePrint("\n>>> Cleaning DataFrame ...")
    if filter_by_ED:
        from matrix_preprocessing_filterData_module import filterBy_randomBC_EditDistance
        any_df = filterBy_randomBC_EditDistance(any_df, inside_ShS=inside_ShS, ED_rule=ED_treshold)
    if False:
        any_df = "HERE NEW FILTERING METHODS"
    verbosePrint(">>> Done!")
#######################################################################################################################################################

### Compute matrixes ##########################################################################################################
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
###############################################################################################################################

### Output #############################################################################################################################################################
verbosePrint("\n>>> Export matrixes ...")
import os
seqCount_matrix_outPath, ShsCount_matrix_outPath, barcodeCount_matrix_outPath, cellCount_matrix_outPath, fragmentEstimate_matrix_outPath  = None, None, None, None, None
# seqCount matrix
seqCount_matrix_outPath = writeMatrix(seqCount_matrix, os.path.join(matrix_files_outdir, "seqCount_matrix.tsv"), matrix_files_delimiter)
# ShsCount matrix
ShsCount_matrix_outPath = writeMatrix(ShsCount_matrix, os.path.join(matrix_files_outdir, "ShsCount_matrix.tsv"), matrix_files_delimiter)
# barcodeCount matrix
barcodeCount_matrix_outPath = writeMatrix(barcodeCount_matrix, os.path.join(matrix_files_outdir, "barcodeCount_matrix.tsv"), matrix_files_delimiter)
# cellCount matrix
cellCount_matrix_outPath = writeMatrix(cellCount_matrix, os.path.join(matrix_files_outdir, "cellCount_matrix.tsv"), matrix_files_delimiter)
# fragmentEstimate matrix
fragmentEstimate_matrix_outPath = writeMatrix(fragmentEstimate_matrix, os.path.join(matrix_files_outdir, "fragmentEstimate_matrix.tsv"), matrix_files_delimiter)
verbosePrint(">>> Matrix Files Created!")
########################################################################################################################################################################

verbosePrint("\n[END]\n")

#+++++++++++++++++++++++++ END CODE ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

