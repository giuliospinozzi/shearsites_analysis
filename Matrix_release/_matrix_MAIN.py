#!/usr/bin/python
# -*- coding: utf-8 -*-


"""
Created on Thu Oct 13 14:29:40 2016

@author: Stefano Brasca
"""

__author__ = "Stefano Brasca"
__copyright__ = "SR-TIGET"
__credits__ = ["Stefano Brasca", "Andrea Calabria", "Giulio Spinozzi", "Adriano De Marino"]
__version__ = "1.0"
__maintainers__ = "Stefano Brasca, Adriano De Marino"
__email__ = "brasca.stefano@hsr.it, demarino.adriano@hsr.it"
__status__ = "Testing"


#++++++++++++++ Requested Package(s) Import +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

### Check requirements, import configuration and PARSE ARGS
import matrix_configure_module

### Import Basic Functions
humanSorted = matrix_configure_module.humanSorted
verbosePrint = matrix_configure_module.verbosePrint
from matrix_preprocessing_dataSources_module import getLaunchPathDict
from matrix_preprocessing_dataLoading_module import loadData
from matrix_processing_computeMatrixes_module import buildSeqCountMatrix, buildShsCountMatrix, buildBarcodeCountMatrix, buildCellCountMatrix, buildFragmentEstimateMatrix, totalMatrix, collisionMatrix
from matrix_output_module import buildOutputPath, writeMatrix


#++++++++++++++++++++++ Global Vars from configure_module +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

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
matrixesTotal_subfolder = matrix_configure_module.matrixesTotal_subfolder
matrixesCollision_subsubfolder = matrix_configure_module.matrixesCollision_subsubfolder
matrix_files_delimiter = matrix_configure_module.matrix_files_delimiter

### Misc configs
dataset_ID = matrix_configure_module.dataset_ID

sc = matrix_configure_module.searchContamination


#++++++++++++++++++++++++ CODE +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

verbosePrint("\n[START]")

### Print Config ##########################################################################################
verbosePrint("\n>>> Configuring ...")
verbosePrint("    DATA:")
verbosePrint("    * dataset_tuple_list: {x}".format(x=str(dataset_tuple_list)))
verbosePrint("    * drop_headers: {x}".format(x=str(drop_headers)))
verbosePrint("    * compression: {x}".format(x=str(compression)))
verbosePrint("    * searchContamination: {x}".format(x=str(sc)))
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
verbosePrint("    * dataset_ID: {x}".format(x=str(dataset_ID)))
verbosePrint("    * common_output_ground_dir: {x}".format(x=str(common_output_ground_dir)))
verbosePrint("    * matrixesTotal_subfolder: {x}".format(x=str(matrixesTotal_subfolder)))
verbosePrint("    * matrixesCollision_subsubfolder: {x}".format(x=str(matrixesCollision_subsubfolder)))
if matrix_files_delimiter == '\t':
    verbosePrint(r'''    * matrix_files_delimiter: \t''')
else:
    verbosePrint("    * matrix_files_delimiter: {x}".format(x=str(matrix_files_delimiter)))
verbosePrint(">>> Done!")
###########################################################################################################

### Check (or create) matrix_files_outdir ############################################################################
verbosePrint("\n>>> Setting up matrixes outdir ...")
matrix_files_outdir = buildOutputPath(common_output_ground_dir)
verbosePrint(">>> Done!")
verbosePrint("\n>>> Setting up matrixesTotal subfolder ...")
matrixesTotal_outdir = buildOutputPath(matrix_files_outdir, matrixesTotal_subfolder)
verbosePrint("\n>>> Setting up matrixesCollision sub-subfolder ...")
matrixesCollision_outdir = buildOutputPath(matrix_files_outdir, matrixesTotal_subfolder, matrixesCollision_subsubfolder)
verbosePrint(">>> Done!")
######################################################################################################################

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

### Compute matrixes without barcode data ##################################################################################################################################
seqCount_matrix, total_seqCount_matrix, ShsCount_matrix, total_ShsCount_matrix, fragmentEstimate_matrix, total_fragmentEstimate_matrix = None, None, None, None, None, None
verbosePrint("\n>>> Computing matrixes:")
verbosePrint("> seqCount matrix ...")
seqCount_matrix = buildSeqCountMatrix(any_df)
total_seqCount_matrix = totalMatrix(seqCount_matrix, column_label=dataset_ID)
collision_seqCount_matrix = collisionMatrix(seqCount_matrix, dataset_ID)
verbosePrint("> ShsCount matrix ...")
ShsCount_matrix = buildShsCountMatrix(any_df)
total_ShsCount_matrix = totalMatrix(ShsCount_matrix, column_label=dataset_ID)
collision_ShsCount_matrix = collisionMatrix(ShsCount_matrix, dataset_ID)
verbosePrint("> fragmentEstimate matrix ...")
fragmentEstimate_matrix = buildFragmentEstimateMatrix(any_df)
total_fragmentEstimate_matrix = totalMatrix(fragmentEstimate_matrix, column_label=dataset_ID)
collision_fragmentEstimate_matrix = collisionMatrix(fragmentEstimate_matrix, dataset_ID)
verbosePrint(">>> Done!")
############################################################################################################################################################################

### Filter Data #######################################################################################################################################
if filter_data:
    verbosePrint("\n>>> Cleaning DataFrame ...")
    if filter_by_ED:
        from matrix_preprocessing_filterData_module import filterBy_randomBC_EditDistance
        any_df = filterBy_randomBC_EditDistance(any_df, inside_ShS=inside_ShS, ED_rule=ED_treshold)
    if False:
        any_df = "HERE NEW FILTERING METHODS"
    verbosePrint(">>> Done!")

#### EDIT BY Adriano and Giulio
if sc:
    import pandas as pd
    import os
    verbosePrint(">>> writing contamination file...")
    any_df.to_csv(os.path.join(matrix_files_outdir,"{dataset_ID}_contamination.csv.gz".format(dataset_ID=str(dataset_ID))), sep='\t', index=False, compression='gzip') #write contamination.csv, this file will be analyze with new function in new program called SKIP.
    verbosePrint(">>> Done!")
#######################################################################################################################################################

### Compute matrixes with barcode data ########################################################################################
barcodeCount_matrix, total_barcodeCount_matrix, cellCount_matrix, total_cellCount_matrix = None, None, None, None
verbosePrint("\n>>> Computing matrixes:")
verbosePrint("> barcodeCount matrix ...")
barcodeCount_matrix = buildBarcodeCountMatrix(any_df)
total_barcodeCount_matrix = totalMatrix(barcodeCount_matrix, column_label=dataset_ID)
collision_barcodeCount_matrix = collisionMatrix(barcodeCount_matrix, dataset_ID)
verbosePrint("> cellCount matrix ...")
cellCount_matrix = buildCellCountMatrix(any_df)
total_cellCount_matrix = totalMatrix(cellCount_matrix, column_label=dataset_ID)
collision_cellCount_matrix = collisionMatrix(cellCount_matrix, dataset_ID)
verbosePrint(">>> Done!")
###############################################################################################################################

### Output ####################################################################################################################################################################################################################################
verbosePrint("\n>>> Export matrixes ...")
import os
seqCount_matrix_outPath, ShsCount_matrix_outPath, barcodeCount_matrix_outPath, cellCount_matrix_outPath, fragmentEstimate_matrix_outPath  = None, None, None, None, None
total_seqCount_matrix_outPath, total_ShsCount_matrix_outPath, total_fragmentEstimate_matrix_outPath, total_barcodeCount_matrix_outPath, total_cellCount_matrix_outPath = None, None, None, None, None
collision_seqCount_matrix_outPath, collision_ShsCount_matrix_outPath, collision_fragmentEstimate_matrix_outPath, collision_barcodeCount_matrix_outPath, collision_cellCount_matrix_outPath = None, None, None, None, None
# seqCount matrix
seqCount_matrix_outPath = writeMatrix(seqCount_matrix, os.path.join(matrix_files_outdir, "{dataset_ID}_seqCount_matrix.tsv".format(dataset_ID=str(dataset_ID))), matrix_files_delimiter)
total_seqCount_matrix_outPath = writeMatrix(total_seqCount_matrix, os.path.join(matrixesTotal_outdir, "total_{dataset_ID}_seqCount_matrix.tsv".format(dataset_ID=str(dataset_ID))), matrix_files_delimiter)
collision_seqCount_matrix_outPath = writeMatrix(collision_seqCount_matrix, os.path.join(matrixesCollision_outdir, "total_{dataset_ID}_seqCount_matrix.tsv".format(dataset_ID=str(dataset_ID))), matrix_files_delimiter)
# ShsCount matrix
ShsCount_matrix_outPath = writeMatrix(ShsCount_matrix, os.path.join(matrix_files_outdir, "{dataset_ID}_ShsCount_matrix.tsv".format(dataset_ID=str(dataset_ID))), matrix_files_delimiter)
total_ShsCount_matrix_outPath = writeMatrix(total_ShsCount_matrix, os.path.join(matrixesTotal_outdir, "total_{dataset_ID}_ShsCount_matrix.tsv".format(dataset_ID=str(dataset_ID))), matrix_files_delimiter)
collision_ShsCount_matrix_outPath = writeMatrix(collision_ShsCount_matrix, os.path.join(matrixesCollision_outdir, "total_{dataset_ID}_ShsCount_matrix.tsv".format(dataset_ID=str(dataset_ID))), matrix_files_delimiter)
# fragmentEstimate matrix
fragmentEstimate_matrix_outPath = writeMatrix(fragmentEstimate_matrix, os.path.join(matrix_files_outdir, "{dataset_ID}_fragmentEstimate_matrix.tsv".format(dataset_ID=str(dataset_ID))), matrix_files_delimiter)
total_fragmentEstimate_matrix_outPath = writeMatrix(total_fragmentEstimate_matrix, os.path.join(matrixesTotal_outdir, "total_{dataset_ID}_fragmentEstimate_matrix.tsv".format(dataset_ID=str(dataset_ID))), matrix_files_delimiter)
collision_fragmentEstimate_matrix_outPath = writeMatrix(collision_fragmentEstimate_matrix, os.path.join(matrixesCollision_outdir, "total_{dataset_ID}_fragmentEstimate_matrix.tsv".format(dataset_ID=str(dataset_ID))), matrix_files_delimiter)
# barcodeCount matrix
barcodeCount_matrix_outPath = writeMatrix(barcodeCount_matrix, os.path.join(matrix_files_outdir, "{dataset_ID}_barcodeCount_matrix.tsv".format(dataset_ID=str(dataset_ID))), matrix_files_delimiter)
total_barcodeCount_matrix_outPath = writeMatrix(total_barcodeCount_matrix, os.path.join(matrixesTotal_outdir, "total_{dataset_ID}_barcodeCount_matrix.tsv".format(dataset_ID=str(dataset_ID))), matrix_files_delimiter)
collision_barcodeCount_matrix_outPath = writeMatrix(collision_barcodeCount_matrix, os.path.join(matrixesCollision_outdir, "total_{dataset_ID}_barcodeCount_matrix.tsv".format(dataset_ID=str(dataset_ID))), matrix_files_delimiter)
# cellCount matrix
cellCount_matrix_outPath = writeMatrix(cellCount_matrix, os.path.join(matrix_files_outdir, "{dataset_ID}_cellCount_matrix.tsv".format(dataset_ID=str(dataset_ID))), matrix_files_delimiter)
total_cellCount_matrix_outPath = writeMatrix(total_cellCount_matrix, os.path.join(matrixesTotal_outdir, "total_{dataset_ID}_cellCount_matrix.tsv".format(dataset_ID=str(dataset_ID))), matrix_files_delimiter)
collision_cellCount_matrix_outPath = writeMatrix(collision_cellCount_matrix, os.path.join(matrixesCollision_outdir, "total_{dataset_ID}_cellCount_matrix.tsv".format(dataset_ID=str(dataset_ID))), matrix_files_delimiter)

verbosePrint(">>> Matrix Files Created!")
###############################################################################################################################################################################################################################################

verbosePrint("\n[END]\n")

#+++++++++++++++++++++++++ END CODE ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

