#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Created on Thu Nov 22 11:29:23 2017

@author: Adriano De Marino
"""

__author__ = "Adriano De Marino"
__copyright__ = "SR-TIGET"
__credits__ = ["Stefano Brasca", "Andrea Calabria", "Giulio Spinozzi","Adriano De Marino"]
__version__ = "1.0"
__maintainer__ = "Stefano Brasca, Adriano De Marino"
__email__ = "demarino.adriano@hsr.it"
__status__ = "Testing"

#++++++++++++++ Requested Package(s) Import +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

### Check requirements, import configuration and PARSE ARGS
import matrix_configure_module

### Import Basic Functions
humanSorted = matrix_configure_module.humanSorted
verbosePrint = matrix_configure_module.verbosePrint
from matrix_preprocessing_dataSources_module import getLaunchPathDict
from matrix_contamination_module import pivot, contamina, comparison, describe, isSymmetric, skip_diag_strided, associationFILE
from matrix_preprocessing_dataLoading_module import loadData
from matrix_processing_computeMatrixes_module import buildSeqCountMatrix, buildShsCountMatrix, buildBarcodeCountMatrix, buildCellCountMatrix, buildFragmentEstimateMatrix, totalMatrix, collisionMatrix
from matrix_output_module import buildOutputPath, writeMatrix, writeMatrixContamination, memory_usage

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
contamination_subfolder = matrix_configure_module.contamination_subfolder
matrixesCollision_subsubfolder = matrix_configure_module.matrixesCollision_subsubfolder
matrix_files_delimiter = matrix_configure_module.matrix_files_delimiter

### Misc configs
dataset_ID = matrix_configure_module.dataset_ID

### --searchContamination option
searchContamination = matrix_configure_module.searchContamination

### --searchContamination option
projectID = matrix_configure_module.projectID

### --metadata option
metadata = matrix_configure_module.metadata

### --onlyContamination option
onlyContamination = matrix_configure_module.onlyContamination

#++++++++++++++++++++++++ CODE +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

verbosePrint("\n[START]")

### Print Config ##########################################################################################
verbosePrint("\n>>> Configuring ...")
verbosePrint("    DATA:")
verbosePrint("    * dataset_tuple_list: {x}".format(x=str(dataset_tuple_list)))
verbosePrint("    * drop_headers: {x}".format(x=str(drop_headers)))
verbosePrint("    * compression: {x}".format(x=str(compression)))
verbosePrint("    * searchContamination: [ {x} ]".format(x=str(searchContamination)))
verbosePrint("    * projectID: [ {x} ]".format(x=str(projectID)))
verbosePrint("    * metadata_file: [ {x} ]".format(x=str(metadata)))
verbosePrint("    * onlyContamination: [ {x} ]".format(x=str(onlyContamination)))
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
if searchContamination:
    verbosePrint("    * contamination_subfolder: {x}".format(x=str(contamination_subfolder)))
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
if searchContamination:
    contamination_outdir = buildOutputPath(matrix_files_outdir, contamination_subfolder)
verbosePrint(">>> Done!")
######################################################################################################################

### Load Data ################################################
launch_path_dict = getLaunchPathDict(dataset_tuple_list)
any_df = loadData(launch_path_dict, drop_headers, compression)
verbosePrint("> Memory usage = {x} Mb".format(x=str(memory_usage())))
##############################################################

### Compute ISs ################################################################################################################################################################
if do_ISs:
    verbosePrint("\n>>> Computing ISs ...")
    from matrix_processing_ISsMethods_module import compute_ISs
    any_df = compute_ISs(any_df, ensembles_per_sample=ensembles_per_sample, ensembles_max_dist=ensembles_max_dist, ensembles_max_span=ensembles_max_span, ISs_method=ISs_method)
    verbosePrint(">>> Done!")
verbosePrint("> Memory usage = {x} Mb".format(x=str(memory_usage())))
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
verbosePrint("> Memory usage = {x} Mb".format(x=str(memory_usage())))
#######################################################################################################################################################


### Output ####################################################################################################################################################################################################################################
verbosePrint("\n>>> Export Contamination matrixes ...")
import os
# Contamination matrix
if searchContamination:
    
    matrix_UMI = pivot(any_df)
    metadata = associationFILE(metadata, projectID)
    df = matrix_UMI.groupby(['chr', 'locus', 'strand', 'shearsite', 'randomBC']).sum()
    matrix = comparison(df)
    summary = describe(matrix,df)
    contamina = contamina(df,metadata)
    
    matrix_UMI_outPath = writeMatrixContamination(matrix_UMI, os.path.join(contamination_outdir, "{dataset_ID}_UMI_matrix.tsv".format(dataset_ID=str(dataset_ID))), matrix_files_delimiter)
    contamina_outPath = writeMatrixContamination(contamina, os.path.join(contamination_outdir, "{dataset_ID}_contamina.tsv".format(dataset_ID=str(dataset_ID))), matrix_files_delimiter)
    summary_outPath = writeMatrixContamination(summary, os.path.join(contamination_outdir, "{dataset_ID}_summary.tsv".format(dataset_ID=str(dataset_ID))), matrix_files_delimiter)


verbosePrint(">>> Matrix Files Created!")
###############################################################################################################################################################################################################################################

verbosePrint("> Memory usage = {x} Mb".format(x=str(memory_usage())))
verbosePrint("\n[END]\n")

exit()
#+++++++++++++++++++++++++ END CODE ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#


