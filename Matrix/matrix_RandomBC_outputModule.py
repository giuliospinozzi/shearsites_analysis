# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 13:46:57 2016

@author: stefano
"""

#++++++++++++++ Requested Package(s) Import +++++++++++++++#
import os, sys
import collections
import matrix_RandomBC_globModule


#++++++++++++++++++++++ Global Vars +++++++++++++++++++++++#
verbose = matrix_RandomBC_globModule.verbose
use_fields = matrix_RandomBC_globModule.use_fields
concat = matrix_RandomBC_globModule.concat


#++++++++++++++++++++++ Global Funcs ++++++++++++++++++++++#
verbosePrint = matrix_RandomBC_globModule.verbosePrint
# humanSorted = matrix_RandomBC_globModule.humanSorted


#+++++++++++++++++++++++++++++++++++++++ FUNCTIONS +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

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


def relabelling(df, asso_dict, use_fields=use_fields, concat=concat, inplace=False):
    # NOTE: - df col names must be the keys of asso_dict!
    #       - use_fields can be an int or a sequence of ints
    
    def buildRelabellingDict(asso_dict, use_fields, concat):
        relabellingDict = {}
        for k, d in asso_dict.items():
            fields = []
            [fields.append(d[n]) for n in use_fields]
            new_label = concat.join(fields)
            relabellingDict[k] = new_label
        return relabellingDict
    
    if isinstance(use_fields, collections.Iterable):
        use_fields = tuple(use_fields)
    else:
        use_fields = tuple((use_fields,))
    
    if inplace is True:
        # copy : boolean, default True
        df.rename(columns=buildRelabellingDict(asso_dict, use_fields, concat), inplace=True)
    else:
        # new object
        return df.rename(columns=buildRelabellingDict(asso_dict, use_fields, concat), inplace=False)
                

def writeMatrix(df, complete_path, out_files_delimiter, metadata=None, verbose=verbose):
    # metadata takes None or asso_dict for relabelling.
    # relabelling function has use_fields and concat fixed as global vars!
    if metadata is None:
        df.to_csv(path_or_buf=complete_path, sep=out_files_delimiter, index_label= 'IS_genomicID', encoding='utf-8')
    else:
        df_relabelled = relabelling(df, metadata)  # see default kwargs
        df_relabelled.to_csv(path_or_buf=complete_path, sep=out_files_delimiter, index_label= 'IS_genomicID', encoding='utf-8')
    verbosePrint(">>> Created file '{complete_path}'".format(complete_path=str(complete_path)))
    return complete_path

#++++++++++++++++++++++++++++++++++++++ MAIN and TEST ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

if __name__ == "__main__":
    pass
  