# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 13:46:57 2016

@author: stefano
"""

#++++++++++++++ Requested Package(s) Import +++++++++++++++#
import os, sys
import matrix_configure_module


#++++++++++++++++++++++ Global Vars +++++++++++++++++++++++#
#verbose = matrix_configure_module.verbose


#++++++++++++++++++++++ Global Funcs ++++++++++++++++++++++#
verbosePrint = matrix_configure_module.verbosePrint
# humanSorted = matrix_RandomBC_globModule.humanSorted


#+++++++++++++++++++++++++++++++++++++++ FUNCTIONS +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#


def buildOutputPath(ground_dir, *subfolders):
    '''
    Purpose: check (and create) an absolute folder path where WRITE FILES
    
    IN: 
    ground_dir - absolute folder path (checked)
    *subfolders - optional, further subfolders to concatenate (checked every time)
    
    OUT:
    OUTDIR - absolute folder path where write, permissions are checked
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
            if sf != '':
                verbosePrint("> subfolder found: {OUTDIR}".format(OUTDIR=str(OUTDIR)))
        # Check: write permissions
        if not os.access(OUTDIR, os.W_OK):
            print "\n[ERROR] You don't have write permissions in OUTDIR='{OUTDIR}'".format(OUTDIR=str(OUTDIR))
            sys.exit("\n[QUIT]\n")
    # Return OUTDIR
    verbosePrint("> OUTDIR: {OUTDIR}".format(OUTDIR=str(OUTDIR)))
    verbosePrint("> write permission: OK.")
    return OUTDIR


def writeMatrix(df, complete_path, out_files_delimiter):
    # df may be None
    if df is None:
        verbosePrint("[Warning] No data available to create '{complete_path}'".format(complete_path=str(complete_path)))
        return None
    df.to_csv(path_or_buf=complete_path, sep=out_files_delimiter, index_label= 'IS_genomicID', encoding='utf-8')
    verbosePrint("> File created: '{complete_path}'".format(complete_path=str(complete_path)))
    return complete_path

def writeMatrixContamination(df, complete_path, out_files_delimiter):
    # df may be None
    if df is None:
        verbosePrint("[Warning] No data available to create '{complete_path}'".format(complete_path=str(complete_path)))
        return None
    df.to_csv(path_or_buf=complete_path, sep=out_files_delimiter, index=True)
    verbosePrint("> File created: '{complete_path}'".format(complete_path=str(complete_path)))
    return complete_path

def memory_usage():
    import resource
    import platform
    memory = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    system = platform.system()

    if system == "Linux" or system == "linux2":
        usage_mb = memory/1024
        return "{:03.2f} MB".format(usage_mb)
    elif system == "Darwin" or system == "darwin":
        usage_mb = memory/1024 ** 2
        return "{:03.2f} MB".format(usage_mb)
