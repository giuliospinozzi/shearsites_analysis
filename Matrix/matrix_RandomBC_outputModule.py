# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 13:46:57 2016

@author: stefano
"""

#++++++++++++++ Requested Package(s) Import +++++++++++++++#
import os, sys
import re

#++++++++++++++++++++++ Global Vars +++++++++++++++++++++++#
import matrix_RandomBC_globModule
verbose = matrix_RandomBC_globModule.verbose

#+++++++++++++++++++++++++++++++++++++++ FUNCTIONS +++++++++++++++++++++++++++++++++++++++#

def verbosePrint(x, verbose=verbose):
    if verbose:
        print x

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


def buildOutputPath(ground_dir, DISEASE, PATIENT, POOL, outfolder):
    # Check ground_dir
    try:
        ground_dir = os.path.normpath(ground_dir)
    except Exception, err_message:
        print "\n[ERROR] ground_dir='{ground_dir}' is not formatted as valid path!".format(ground_dir=str(ground_dir))
        print "os.path.normpath returned: ", err_message
        sys.exit("\n[QUIT]\n")
    # Generate OUTDIR path
    BASEDIR = os.path.normpath(os.path.join(ground_dir, DISEASE, PATIENT))
    OUTDIR = os.path.normpath(os.path.join(BASEDIR, "quantification", POOL, "RandomBC", outfolder))
    # Create and/or check
    if os.path.isfile(OUTDIR):
        print "\n[ERROR] OUTDIR must be a folder path, not a file! Check buildOutputPath parameters; OUTDIR: '{OUTDIR}'".format(OUTDIR=str(OUTDIR))
        sys.exit("\n[QUIT]\n")
    if not os.path.exists(OUTDIR):
        try:
            os.makedirs(OUTDIR)
        except Exception, err_message:
            print "\n[ERROR] OUTDIR is not a valid path or you don't have sufficient privileges to create it! OUTDIR: '{OUTDIR}'".format(OUTDIR=str(OUTDIR))
            print "os.makedirs returned: ", err_message
            sys.exit("\n[QUIT]\n")
        verbosePrint("> OUTDIR created: {OUTDIR}".format(OUTDIR=str(OUTDIR)))
    else:
        verbosePrint("> OUTDIR found: {OUTDIR}".format(OUTDIR=str(OUTDIR)))
    if not os.access(OUTDIR, os.W_OK):
        print "\n[ERROR] Can't write anything in OUTDIR: '{OUTDIR}'".format(OUTDIR=str(OUTDIR))
        sys.exit("\n[QUIT]\n")
    return OUTDIR

   
def writeMatrix(df, complete_path, out_files_delimiter, metadata=None, verbose=verbose):
    # metadata is thought to get asso_dict as input
    df.to_csv(path_or_buf=complete_path, sep=out_files_delimiter, index_label= 'IS_genomicID', encoding='utf-8')
    verbosePrint("> Created file '{complete_path}'".format(complete_path=str(complete_path)))
    return complete_path

#++++++++++++++++++++++++++++++++++++++ MAIN and TEST ++++++++++++++++++++++++++++++++++++++#

if __name__ == "__main__":
    pass
  