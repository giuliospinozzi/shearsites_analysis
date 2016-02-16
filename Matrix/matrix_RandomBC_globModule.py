# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 09:49:24 2016

@author: stefano
"""


#++++++++++++++++++++++++++++++++++ Global Vars +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

### SHARED VARS ####################################################################################################################################

## Screen print ##
verbose = True

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

## Output - outputModule ##
outfolder = "Matrixes"  #ground_dir, DISEASE, PATIENT, POOL as for Data -> out path = in path + outfolder
out_files_delimiter = '\t'
# Relabel columns (BARCODES <-> master-keys of asso data dict in this implementation <-> sample column in any_df)
# with a concatenation of related attributes (fields/columns of AssoFile <-> sub-keys of asso data dict sub-dicts)
relabelling = True
use_fields = 6  # can be an int or a sequence of ints
concat = "_"  # char to concatenate fields in use_fields

####################################################################################################################################################


### DIAGNOSTIC VARS ################################################################################################################################

## Export MATRIXES ##
# DO
export_matrixes = True

## Export CEM data ##
# DO
export_cem_data = True
cem_data_outfolder = "CEMdata"  #ground_dir, DISEASE, PATIENT, POOL as for Data -> CEM data out path = in path + cem_data_outfolder
# build out file name
data_id = ""  # something related to data_files_name_filter
out_filename = "CEMdata" + "_" + DISEASE + "_" + POOL
if data_id: out_filename += "_" + data_id 
out_filename += ".tsv"

## Export DIAGNOSTICS ##
# DO
export_diagnostics = False
diagnostic_outfolder = "Diagnostics"  #ground_dir, DISEASE, PATIENT, POOL as for Data -> diagnostic out path = in path + diagnostic_outfolder
# Data selection
specific_samples = True  # 'True' here requires explicit lists below
dilution_to_process = ['L', 'M', 'N']  # Whatever, if specific_samples is False
condition_to_process = ['LMv2-II']  # Whatever, if specific_samples is False
# Task to perform
checkNucleotidesBalancing = True  # Stacked-bar plot 
FragmentLengthDistribution = True  # Fragment length histogram, fitted with soncLength and KDE
checkShearSitesOccurrency = True  # Occurrency bar plot
checkRandomBCoccurrency = True  # Occurrency line plot
checkEditDistance_diagonal = True   # Edit Distance occurrency histogram (within shearsites)
checkEditDistance_extensive = False   # Edit Distance occurrency histogram  (all-VS-all)
plot_heatmap = True  # evaluated only if checkEditDistance_XXX is True  # Edit distance heatmap
limit_heatmap_plot = (True, 600)  # or (False, whatever); syntax: (Do?, max number of rows-cols allowed)
                                  # evaluated only if plot_heatmap is True
####################################################################################################################################################


#++++++++++++++++++++++++++++++++++ Global Funcs +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

def verbosePrint(x, verbose=verbose):
    if verbose:
        print x

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
