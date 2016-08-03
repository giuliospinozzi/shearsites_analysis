# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 09:49:24 2016

@author: stefano
"""


#++++++++++++++++++++++++++++++++++ Global Vars +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

## Screen print ##
verbose = True
print_time = True # evaluated only if verbose is True

## Association File - matrix_preprocessing_assoFile_module ##
asso_folder = "/opt/applications/scripts/isatk/elements/association"
asso_file_name = "asso.assayvalidation.lane2.tsv"
asso_delimiter = '\t'

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

### Filter Data - matrix_preprocessing_filterData_module ##
filter_data = True  # allow filtering as stated below
# filterBy_randomBC_EditDistance: vars evaluated (and tasks executed) only if filter_data is True
filter_by_ED = True
inside_ShS = True
ED_treshold = 3  # int (1->11); can be even a callable (e.g. ED_rule=func such that func(arg)-> int. See matrix_preprocessing_filterData_module)

## COMMON OUTPUT GROUND DIR ##
common_output_ground_dir = "/storage/d3/tmp/stefano/test_Matrix"  # abs path used as 'ground dir' for subfolder tree.

## Matrix output Files - outputModule ##
matrix_outfolder = "Matrixes"  # subfolder of common_output_ground_dir where write matrixes
matrix_files_delimiter = '\t'
# Relabel columns (BARCODES <-> master-keys of asso data dict in this implementation <-> sample column in any_df)
# with a concatenation of related attributes (fields/columns of AssoFile <-> sub-keys of asso data dict sub-dicts)
relabelling = True
use_fields = 6  # can be an int or a sequence of ints
concat = "_"  # char to concatenate fields in use_fields (if more than one)



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

# OBSOLETE
#==============================================================================
# import collections   
# def flattenDict(d, parent_key='', sep='@@'):
#     '''
#     - sep is used only in computation, so a 'robust' separator is advised
#     - parent_key can be used to add a fixed string on the key begins while flattening
#     '''
#     items = []
#     for k, v in d.items():
#         new_key = parent_key + sep + str(k) if parent_key else k
#         if isinstance(v, collections.MutableMapping):
#             items.extend(flattenDict(v, new_key, sep=sep).items())
#         else:
#             items.append((tuple(new_key.split(sep)), v))
#     return dict(items)
#==============================================================================
# OBSOLETE
#==============================================================================
# import sys
# def verbosePrint(x, verbose=verbose):
#     if verbose:
#         print x
#         sys.stdout.flush()
#==============================================================================


#++++++++++++++++++++++++++++++++++++++ MAIN and TEST +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

if __name__ == "__main__":
    pass
