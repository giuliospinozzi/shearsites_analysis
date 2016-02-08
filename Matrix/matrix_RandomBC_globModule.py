# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 09:49:24 2016

@author: stefano
"""


#++++++++++++++++++++++++++++++++++ Global Vars ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

### SHARED VARS ################################################################

# Screen print
verbose = True

# Association File - assoModule
asso_folder = "/opt/applications/scripts/isatk/elements/association"
asso_file_name = "asso.assayvalidation.lane1.tsv"
asso_delimiter = '\t'

# Data - dataModule
ground_dir = "/opt/NGS/results"
DISEASE = "AssayValidation"
PATIENT = "CEMJY"
POOL = "LANE_1"
data_files_delimiter = '\t'
data_files_name_filter = ".ISfixed.LENGTHfixed.randomBC.tsv"  # BE CAREFUL HERE!

# Output - outputModule
#ground_dir, DISEASE, PATIENT, POOL as for Data
outfolder = "Matrixes"
out_files_delimiter = '\t'
relabelling = True
use_fields = 6  # can be an int or a sequence of ints
concat = "_"  # char to concatenate fields

################################################################################


### DIAGNOSTIC VARS ############################################################

# Export MATRIXES
export_matrixes = False

# Export CEM data
# config
export_cem_data = False
data_id = ""  # something related to data_files_name_filter
# build out file name structure
out_filename = "CEMdata" + "_" + DISEASE + "_" + POOL
if data_id: out_filename += "_" + data_id 
out_filename += ".tsv"

# Export DIAGNOSTICS
export_diagnostics = True

################################################################################


#++++++++++++++++++++++++++++++++++ Global Funcs ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

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


#++++++++++++++++++++++++++++++++++++++ MAIN and TEST +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

if __name__ == "__main__":
    pass
