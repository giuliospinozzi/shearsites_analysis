# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 09:49:24 2016

@author: stefano
"""


#++++++++++++++++++++++ Global Vars +++++++++++++++++++++++#

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
data_files_name_filter = ".randomBC.tsv"  # always valid despite corrections applied

# Output - outputModule
#ground_dir, DISEASE, PATIENT, POOL as for Data
outfolder = "Matrixes"
out_files_delimiter = '\t'
relabelling = True
use_fields = 6  # can be an int or a sequence of ints
concat = "_"  # char to concatenate fields

# Export CEM
export_cem = True

#++++++++++++++++++++++ Global Funcs +++++++++++++++++++++++#

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

#++++++++++++++++++++++++++++++++++++++ MAIN and TEST ++++++++++++++++++++++++++++++++++++++#

if __name__ == "__main__":
    pass
