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


def buildOutputPath(ground_dir, *subfolders):
    '''
    Purpose: check (and create) an absolute folder path where WRITE FILES
    
    IN: 
    ground_dir - absolute folder path (checked)
    *subfolders - optional, further subfolders to concatenate (checked every time)
    
    OUT:
    OUTDIR - absolute folder path where write permissions are checked
    '''
    # Final path variable
    OUTDIR = None  # <- os.path.normpath(ground_dir) + join *subfolders
    # Try OUTDIR=normpath(OUTDIR) and formal check
    try:
        OUTDIR = os.path.normpath(ground_dir)
    except Exception, err_message:
        print "\n[ERROR] ground_dir='{ground_dir}' is not formatted as valid path!".format(ground_dir=str(ground_dir))
        print "os.path.normpath returned: ", err_message
        sys.exit("\n[QUIT]\n")
    if os.path.isfile(OUTDIR):
        print "\n[ERROR] ground_dir='{OUTDIR}' must be a folder path, not a file!".format(OUTDIR=str(OUTDIR))
        sys.exit("\n[QUIT]\n")
    if not os.path.isabs(OUTDIR):
        print "\n[ERROR] ground_dir='{OUTDIR}' must be an absolute path!".format(OUTDIR=str(OUTDIR))
        sys.exit("\n[QUIT]\n")
    # Check ground_dir: create it if not exists
    if not os.path.exists(OUTDIR):
        try:
            os.makedirs(OUTDIR)
        except Exception, err_message:
            print "\n[ERROR] ground_dir='{OUTDIR}' is not a valid path or you don't have permissions to create it!".format(OUTDIR=str(OUTDIR))
            print "os.makedirs returned: ", err_message
            sys.exit("\n[QUIT]\n")
        verbosePrint("> ground_dir created: {OUTDIR}".format(OUTDIR=str(OUTDIR)))
    else:
        verbosePrint("> ground_dir found: {OUTDIR}".format(OUTDIR=str(OUTDIR)))
    # Check ground_dir: write permissions if required
    if not subfolders:
        if not os.access(OUTDIR, os.W_OK):
            print "\n[ERROR] You don't have write permissions in ground_dir='{OUTDIR}'".format(OUTDIR=str(OUTDIR))
            sys.exit("\n[QUIT]\n")
    # Loop over *subfolders
    for sf in subfolders:
        # Try OUTDIR=join(OUTDIR, sf) and formal check
        try:
            OUTDIR = os.path.normpath(os.path.join(OUTDIR, sf))
        except Exception, err_message:
            print "\n[ERROR] Cannot join OUTDIR='{OUTDIR}' and SUBFOLDER='{SUBFOLDER}' as a valid path!".format(OUTDIR=str(OUTDIR), SUBFOLDER=str(sf))
            print "os.path.normpath(os.path.join(...)) returned: ", err_message
            sys.exit("\n[QUIT]\n")
        if os.path.isfile(OUTDIR):
            print "\n[ERROR] *subfolders args must be folder(s), not file(s)! Input subfolders: {subfolders}.".format(subfolders=str(subfolders))
            sys.exit("\n[QUIT]\n")
        # Check: create OUTDIR if not exists
        if not os.path.exists(OUTDIR):
            try:
                os.makedirs(OUTDIR)
            except Exception, err_message:
                print "\n[ERROR] OUTDIR='{OUTDIR}' + SUBFOLDER='{SUBFOLDER}' is not a valid path or you don't have permissions to create it!".format(OUTDIR=str(OUTDIR), SUBFOLDER=str(sf))
                print "os.makedirs returned: ", err_message
                sys.exit("\n[QUIT]\n")
            verbosePrint("> subfolder created: {OUTDIR}".format(OUTDIR=str(OUTDIR)))
        else:
            verbosePrint("> subfolder found: {OUTDIR}".format(OUTDIR=str(OUTDIR)))
        # Check: write permissions
        if not os.access(OUTDIR, os.W_OK):
            print "\n[ERROR] You don't have write permissions in OUTDIR='{OUTDIR}'".format(OUTDIR=str(OUTDIR))
            sys.exit("\n[QUIT]\n")
    # Return OUTDIR
    verbosePrint(">>> OUTDIR: {OUTDIR}".format(OUTDIR=str(OUTDIR)))
    verbosePrint(">>> OUTDIR WRITE PERMISSIONS: [OK]")
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
    verbosePrint(">>> File created: '{complete_path}'".format(complete_path=str(complete_path)))
    return complete_path

#++++++++++++++++++++++++++++++++++++++ MAIN and TEST ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

if __name__ == "__main__":
    pass
  