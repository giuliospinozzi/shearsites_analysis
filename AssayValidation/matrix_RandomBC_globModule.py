# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 09:49:24 2016

@author: stefano
"""


#++++++++++++++++++++++++++++++++++ Global Vars +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

### SHARED VARS ####################################################################################################################################

## Screen print ##
verbose = True
print_time = True # evaluated only if verbose is True

## Association File - assoModule ##
asso_folder = "/opt/applications/scripts/isatk/elements/association"
asso_file_name = "asso.assayvalidation.lane2.tsv"
asso_delimiter = '\t'

## Data - dataModule ##
ground_dir = "/opt/NGS/results"
DISEASE = "AssayValidation"
PATIENT = "CEMJY"
POOL = "LANE_2"
data_files_delimiter = '\t'
data_files_name_filter = ".ISfixed.LENGTHfixed.randomBC.tsv"  # BE CAREFUL HERE!

## COMMON OUTPUT GROUND DIR ##
common_output_ground_dir = "/storage/d3/tmp/stefano/randomBC_quantification_filterByHeaders/{POOL}".format(POOL=POOL)  # PUT HERE AN ABS PATH

## Matrix output Files - outputModule ##
matrix_outfolder = "Matrixes"
matrix_files_delimiter = '\t'
# Relabel columns (BARCODES <-> master-keys of asso data dict in this implementation <-> sample column in any_df)
# with a concatenation of related attributes (fields/columns of AssoFile <-> sub-keys of asso data dict sub-dicts)
relabelling = True
use_fields = 6  # can be an int or a sequence of ints
concat = "_"  # char to concatenate fields in use_fields

####################################################################################################################################################


### FILTERING CONFIG #############################################
filter_data = False
byHeaders = False  # evaluated only if filter_data is True
bySC = False  # evaluated only if filter_data is True
##################################################################


### DIAGNOSTIC VARS ################################################################################################################################

## Export MATRIXES ##
# DO
export_matrixes = True  # path = common_output_ground_dir+matrix_outfolder

## Export CEM data ##
# DO
export_cem_data = True
cem_data_outfolder = "CEMdata"  # path = common_output_ground_dir+cem_data_outfolder (build in CEMmodule)
# build out file name
data_id = ""  # something related to data_files_name_filter
cem_data_outfile_name = "CEMdata" + "_" + DISEASE + "_" + POOL
if data_id: cem_data_outfile_name += "_" + data_id 
cem_data_outfile_name += ".tsv"

## Export DIAGNOSTICS ##
# DO
export_diagnostics = True
diagnostic_outfolder = "Diagnostics"  # path = common_output_ground_dir+diagnostic_outfolder
# Data selection
specific_samples = True  # 'True' here requires explicit lists below
condition_to_process = ['LMv2-II']  # evaluated only if specific_samples is True
approach_to_process = ['Block']  # evaluated only if specific_samples is True
dilution_to_process = ['L']  # evaluated only if specific_samples is True  # ['L', 'M', 'N']
# Task granularity
inside_ISs = False
# Task to perform
checkNucleotidesBalancing = False  # Stacked-bar plot 
FragmentLengthDistribution = False  # Fragment length histogram, fitted with soncLength and KDE
checkShearSitesOccurrency = False  # Occurrency bar plot
checkRandomBCoccurrency = False  # Occurrency line plot
checkEditDistance_diagonal = False   # Edit Distance occurrency histogram (within shearsites)
checkEditDistance_extensive = True   # Edit Distance occurrency histogram  (all-VS-all)
plot_heatmap = False  # evaluated only if checkEditDistance_XXX is True  # Edit distance heatmap
limit_heatmap_plot = (True, 800)  # or (False, whatever); syntax: (Do?, max number of rows-cols allowed)
                                  # evaluated only if plot_heatmap is True
plot_heatmap_byChunks = False  # evaluated only if checkEditDistance_XXX is True  # Many Edit distance sub-heatmap
ShS_chunk_size = 11  # evaluated only if plot_heatmap_byChunks is True; should be int>=3, odd.
checkBCcountRatio = False  # Violin plot of randomBC seq-count ratios, divided in classes by edit-distance
####################################################################################################################################################


#++++++++++++++++++++++++++++++++++ Global Funcs +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#==============================================================================
# import sys
# def verbosePrint(x, verbose=verbose):
#     if verbose:
#         print x
#         sys.stdout.flush()
#==============================================================================

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

import collections   
def flattenDict(d, parent_key='', sep='@@'):
    '''
    - sep is used only in computation, so a 'robust' separator is advised
    - parent_key can be used to add a fixed string on the key begins while flattening
    '''
    items = []
    for k, v in d.items():
        new_key = parent_key + sep + str(k) if parent_key else k
        if isinstance(v, collections.MutableMapping):
            items.extend(flattenDict(v, new_key, sep=sep).items())
        else:
            items.append((tuple(new_key.split(sep)), v))
    return dict(items)


#++++++++++++++++++++++++++++++++++++++ MAIN and TEST +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

if __name__ == "__main__":
    pass
