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
out_files_delimiter = "\t"

#++++++++++++++++++++++ Global Funcs +++++++++++++++++++++++#

def verbosePrint(x, verbose=verbose):
    if verbose:
        print x

#++++++++++++++++++++++++++++++++++++++ MAIN and TEST ++++++++++++++++++++++++++++++++++++++#

if __name__ == "__main__":
    pass
