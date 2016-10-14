# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 09:49:24 2016

@author: stefano
"""


#++++++++++++++++++++++++++++++++++ Global Vars +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

### This section can be exploited to launch matrix_MAIN_MODULE directly, without parsing args
### (Console usage or debug purposes)
### In this case, the section 'Override Global Vars by argumentParser_module' below must be commented
### Otherwise all these variables will be overrided!

## Screen print ##
verbose = True
print_time = True # evaluated only if verbose is True

## Data - matrix_preprocessing_dataSources_module, matrix_preprocessing_dataGathering_module, matrix_preprocessing_dataLoading_module ##
dataset_tuple_list = [ ('/opt/NGS/results', 'AssayValidation', 'CEMJY', ['LANE_2']) ]
# dataset_tuple_list is a list of tuple(s) like:
# (abs-path-where-find-data, a-disease, a-patient-of-the-disease, specific-pool-selection?)
#   examples:
#   --> ('/opt/NGS/results', 'AssayValidation', 'CEMJY', False) - all pools available (LANE_1 and LANE_2)
#   or
#   --> ('/opt/NGS/results', 'AssayValidation', 'CEMJY', ['LANE_2']) - only pool(s) in list (just 'LANE_2 in this case)
drop_headers = True  # if False, headers are kept in any_df (as lists, under 'header_list' column)
compression = 'gzip'  # supported compression: None or 'gzip'. For both refactored data and randomBC data.

## ISs computation - matrix_processing_ISsMethods_module ##
do_ISs = True
# ensembles config: vars evaluated (and tasks executed) only if do_ISs is True
ensembles_per_sample = False
ensembles_max_dist = 7
ensembles_max_span = 8
# ISs method: vars evaluated (and tasks executed) only if do_ISs is True
ISs_method = 'classic'  # The only one available up to now. place_on_mode=True default behaviour.

## Filter Data - matrix_preprocessing_filterData_module ##
filter_data = True  # allow filtering as stated below
# filterBy_randomBC_EditDistance: vars evaluated (and tasks executed) only if filter_data is True
filter_by_ED = True
inside_ShS = True
ED_treshold = 3  # int (1->11); can be even a callable (e.g. ED_rule=func such that func(arg)-> int. See matrix_preprocessing_filterData_module)

## COMMON OUTPUT GROUND DIR ##
common_output_ground_dir = "/storage/d3/tmp/stefano/debug"  # abs path used as 'ground dir' for subfolder tree.

## Matrix output Files - outputModule ##
matrix_outfolder = "Matrixes"  # subfolder of common_output_ground_dir where write matrixes
matrix_files_delimiter = '\t'


#++++++++++++++++++++++++++++++++++ Override Global Vars by argumentParser_module ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

### This section overrides variables above according to matrix_argumentParser_module
### and Costants below

## LAUNCH ARGUMENT PARSER ##
import matrix_argumentParser_module
args = matrix_argumentParser_module.args

## Costants
# Note: program behaviour, arg parsing, help, description...
#       the whole program leverage on these 'costants' settings!
#       DO NOT CHANGE ANYTHING.
print_time = True
drop_headers = True
compression = 'gzip'
filter_by_ED = True
ISs_method = 'classic'
matrix_files_delimiter = '\t'

## Screen print
verbose = args.quiet
## Data
dataset_tuple_list = args.dataset_tuple_list
## ISs computation
do_ISs = args.covered_bases
ensembles_per_sample = args.ISs_per_sample
ensembles_max_dist = args.ISs_max_gap
ensembles_max_span = args.ISs_max_span
## Filter Data
filter_data = args.unfiltered
inside_ShS = args.filter_ignoring_shearsites
ED_treshold = args.filter_edit_distance_threshold


#++++++++++++++++++++++++++++++++++ Global Funcs +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

import sys
from time import localtime, strftime
def verbosePrint(x, verbose=verbose, print_time=print_time):
    if verbose:
        if print_time:
            y = str(x)
            nnl = y.count("\n", 0, len(x)/2+1)
            y = y.replace("\n", "", nnl)
            nl_str = "".join(['\n']*nnl)
            print nl_str+"[{time}] ".format(time=strftime("%Y-%m-%d %H:%M:%S", localtime())), y
        else:
            print x
        sys.stdout.flush()

import re
def humanSorted(l):
    def tryint(s):
        try:
            return int(s)
        except:
            return s
    def alphanum_key(s):
        return [ tryint(c) for c in re.split('([0-9]+)', s) ]
    def sort_nicely(l):
        return sorted(l, key=alphanum_key)
    return sort_nicely(l)


#++++++++++++++++++++++++++++++++++ Code +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

### This section prints the whole configuration in order to be sure about the launch
### and avoid uncontrolled overriding!

### MOVED IN MAIN
#verbosePrint("\n>>> Configuring ...")
#verbosePrint("    DATA:")
#verbosePrint("    * dataset_tuple_list: {x}".format(x=str(dataset_tuple_list)))
#verbosePrint("    * drop_headers: {x}".format(x=str(drop_headers)))
#verbosePrint("    * compression: {x}".format(x=str(compression)))
#verbosePrint("    ISs:")
#verbosePrint("    * do_ISs: {x}".format(x=str(do_ISs)))
#if do_ISs:
#    verbosePrint("     * ensembles_per_sample: {x}".format(x=str(ensembles_per_sample)))
#    verbosePrint("     * ensembles_max_dist: {x}".format(x=str(ensembles_max_dist)))
#    verbosePrint("     * ensembles_max_span: {x}".format(x=str(ensembles_max_span)))
#    verbosePrint("     * ISs_method: {x}".format(x=str(ISs_method)))
#verbosePrint("    FILTER:")
#verbosePrint("    * filter_data: {x}".format(x=str(filter_data)))
#if filter_data:
#    verbosePrint("     * filter_by_ED: {x}".format(x=str(filter_by_ED)))
#    verbosePrint("     * inside_ShS: {x}".format(x=str(inside_ShS)))
#    verbosePrint("     * ED_treshold: {x}".format(x=str(ED_treshold)))
#verbosePrint("    OUTPUT:")
#verbosePrint("    * common_output_ground_dir: {x}".format(x=str(common_output_ground_dir)))
#if matrix_outfolder != '':
#    verbosePrint("    * matrix_outfolder: {x}".format(x=str(matrix_outfolder)))
#if matrix_files_delimiter == '\t':
#    verbosePrint(r'''    * matrix_files_delimiter: \t''')
#else:
#    verbosePrint("    * matrix_files_delimiter: {x}".format(x=str(matrix_files_delimiter)))
#verbosePrint(">>> Done!")

