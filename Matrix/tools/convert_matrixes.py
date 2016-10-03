#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 11:06:01 2016

@author: stefano
"""


import numpy as np
import pandas as pd
import sys
import os


__author__ = "Stefano Brasca"
__copyright__ = "SR-TIGET"
__credits__ = ["Stefano Brasca", "Andrea Calabria", "Giulio Spinozzi"]
__version__ = "1.0"
__maintainer__ = "Stefano Brasca"
__email__ = "brasca.stefano@hsr.it"
__status__ = "Testing"


### APPUNTI DI SVILUPPO ##############################################################################################################
# 1) - FATTO !!! l'argomento inpath deve accettare non solo (liste di) file ma anche (liste di) directory, da processare per intero.
#    Deve accorgersi automaticamente se si tratta di directory e, in questo caso, deve
#    occuparsi di ricavare tutti i path dei file da processare e fare override della variabile 'input_paths'
# 2) aggiungere altri arg che ora sono parametri statici: qualcosa è già pronto (es. old_to_new_suffix ...), altri
#    vanno cercati nel codice (es. input_matrixes_encoding ...)
# 3) SISTEMARE HELP e description
######################################################################################################################################


### DEFAULT VALUES AS MODULE #############################################################################################
# Print on screen
verbose = True
print_time = True  # eval if verbose is True
# Output
old_to_new_suffix = '_newFormat.tsv'
new_to_old_suffix = '_oldFormat.tsv'
outdir = None  # write outfile(s) alongside infile(s) [None] or specify a different common folder [complete out dir path]
##########################################################################################################################


if __name__ == '__main__':
    
    description = """
+---+ GENERAL DESCRIPTION +---------------------------------------------------+

Write here a brief description

+---+ DETAILS +---------------------------------------------------------------+

Write here a detailed description

+---+ ABOUT ARGUMENTS: +------------------------------------------------------+ """

    ### DEFAULT VALUES AS PROGRAM ##################################################################################################
    # Print on screen
    verbose = False
    print_time = True  # eval if verbose is True
    # Output
    old_to_new_suffix = '_newFormat.tsv'  # NOT PARSED UP TO NOW
    new_to_old_suffix = '_oldFormat.tsv'  # NOT PARSED UP TO NOW
    outdir = None  # write outfile(s) alongside infile(s) [None] or specify a different common folder [complete out dir path]
    # Command choices
    #command_choices = ['old_to_new', 'new_to_old']
    ################################################################################################################################

    ### ARGS ############################################################################################################################
    import argparse
    class MyParser(argparse.ArgumentParser): 
       def error(self, message):
          print 'error: %s\n' % message
          print ""
          self.print_help()
          sys.exit("\n[QUIT]\n")
    # Parser
    parser = MyParser(description=description, formatter_class=argparse.RawTextHelpFormatter)
    # Shared args
    parser.add_argument("inpath", metavar='INPATHs', nargs='+', help="absolute file path(s) for input matrix(es) (to be extendex with directories' paths ... )")
    parser.add_argument("-od", "--output_directory_path", metavar='OUTDIR', default=outdir, help="write help for OUTDIR ...")
    parser.add_argument("-v", "--verbose", action="store_true", help="verbose execution")
    # SubParsers
    subparsers = parser.add_subparsers(dest='command', metavar='COMMAND', help='write help for COMMAND ...')
    # create the parser for the "old_to_new" command
    parser_old_to_new = subparsers.add_parser('old_to_new', help='write help for "old_to_new" command ...')
    ## Args for "old_to_new"
    #parser_old_to_new.add_argument('--arg_old_to_new', help='arg_old_to_new help')
    # create the parser for the "new_to_old" command
    parser_new_to_old = subparsers.add_parser('new_to_old', help='write help for "new_to_old" command ...')
    ## Args for "new_to_old"
    #parser_new_to_old.add_argument('--arg_new_to_old', help='arg_new_to_old help')
    # Parse Args
    args = parser.parse_args()
    
    #####################################################################################################################################

    ### CHECK AND SET VALUES ############################################################################################################
    
    # TO DO: if a directory (or more) is given, all the files inside will be processed!
    input_paths = args.inpath
    
    if args.verbose != verbose:
        verbose = args.verbose
    if args.output_directory_path != outdir:
        outdir = args.output_directory_path
        
    if args.command == 'old_to_new':
        # here specific args and settings
        # e.g.: old_to_new_suffix
        pass
    elif args.command == 'new_to_old':
        # here specific args and settings
        # e.g.: new_to_old_suffix
        pass
    #####################################################################################################################################


### FUNC #################################################################################################################################################################################

def verbosePrint(x, verbose=verbose, print_time=print_time):
    '''
    Purpose: print immediately if verbose is True, with
             date and time if print_time is True too.
    
    IN: 
    x - anything supporting __str__ method
    
    OUT:
    x on stdout, immediately flushed
    
    NOTE:
    uncommon/complex patterns of "\n" may be not properly
    handled in case of 'print_time'. 
    '''
    if verbose:
        if print_time:
            from time import localtime, strftime
            y = str(x)
            nnl = y.count("\n", 0, len(x)/2+1)
            y = y.replace("\n", "", nnl)
            nl_str = "".join(['\n']*nnl)
            print nl_str+"[{time}] ".format(time=strftime("%Y-%m-%d %H:%M:%S", localtime())), y
        else:
            print x
        sys.stdout.flush()

def check_input_paths (*paths, **kwargs):
    '''
    Purpose: check input matrix file path(s);
             exit and explain why if something is wrong.
    
    IN: 
    *paths - complete path(s) of matrix file(s)
    
    OUT:
    0
    '''
    # set behaviour by kwargs
    quiet=kwargs.get('quiet', False)
    check_isfile=kwargs.get('check_isfile', True)

    if not quiet:
        verbosePrint("\n[CHECK IN PATH(S)]")
    for p in paths:
        # normalize path
        p = str(p)
        try:
            p = os.path.normpath(p)
        except Exception, err_message:
            print "\n[ERROR] input path='{p}' is not formatted as valid path!".format(p=str(p))
            print "os.path.normpath returned: ", err_message
            sys.exit("\n[QUIT]\n")
        p = os.path.expanduser(p)
        p = os.path.expandvars(p)
        # check if path is absolute, read permission and if path points to an existing file
        if not os.path.isabs(p):
            print "\n[ERROR] input path must be an absolute path! Your input path='{p}'".format(p=str(p))
            sys.exit("\n[QUIT]\n")
        if not os.access(os.path.dirname(p), os.R_OK):
            print "\n[ERROR] You have not read permission in folder='{f}'".format(f=str(os.path.dirname(p)))
            sys.exit("\n[QUIT]\n")
        if check_isfile:
            if not os.path.isfile(p):
                print "\n[ERROR] input path must point to an existing file! Your input path='{p}'".format(p=str(p))
                sys.exit("\n[QUIT]\n")
    if not quiet:
        verbosePrint("...OK!")
        verbosePrint("[DONE]")
    return 0

def input_paths_from_dir (*paths):
    '''
    '''
    verbosePrint("\n[GET IN PATH(S)]")
    check_input_paths (*paths, quiet=True, check_isfile=False)
    path_list = []
    for p in paths:
        if os.path.isdir(p):
            path_list += [os.path.normpath(os.path.join(p, f)) for f in os.listdir(p) if (os.path.isfile(os.path.join(p, f)) and not p.startswith('.'))]
        else:
            path_list.append(p)
    verbosePrint("...OK!")
    verbosePrint("* input paths: {path_list}".format(path_list=str(path_list)))
    verbosePrint("[DONE]")
    return path_list

def normalize_input_paths (*paths):
    '''
    Purpose: return normalized paths
    (normpath, expanduser, expandvars)
    after being checked by check_input_paths
    
    IN: 
    *paths - complete path(s) of matrix file(s)
    
    OUT:
    paths_norm - tuple of normalized paths
    '''
    paths_norm = []
    for p in paths:
        p_norm = os.path.normpath(p)
        p_norm = os.path.expanduser(p_norm)
        p_norm = os.path.expandvars(p_norm)
        p_norm = os.path.normpath(p)
        paths_norm.append(p_norm)
    return paths_norm

def import_new_matrixes (*paths):
    '''
    Purpose: load matrix file(s) - new format
    as Pandas DataFrame(s)
    
    IN: 
    *paths - complete path(s) of matrix file(s)
    
    OUT:
    matrixes - list of DataFrame objects
    '''
    # Settings
    input_matrixes_encoding = 'utf-8'
    input_matrixes_sep = '\t'
    # Func
    def parse_matrix (path):
        verbosePrint("> loading {path} ... ".format(path=str(path)))
        return pd.DataFrame.from_csv(path,encoding=input_matrixes_encoding,sep=input_matrixes_sep,parse_dates=False)
    # Call
    verbosePrint("\n[LOAD MATRIX(ES) - New Format]")
    matrixes = [parse_matrix(p) for p in paths]
    verbosePrint("[DONE]")
    return matrixes
    
def import_old_matrixes (*paths):
    '''
    Purpose: load matrix file(s) - old format
    as Pandas DataFrame(s)
    
    IN: 
    *paths - complete path(s) of matrix file(s)
    
    OUT:
    matrixes - list of DataFrame objects
    '''
    # Parse settings
    input_matrixes_encoding = 'utf-8'
    input_matrixes_sep = '\t'
    index_col = (0,1,2)
    # Clean settings
    cols_to_drop = ['all', 'GeneName', 'GeneStrand']
    drop_cols_starting_with = ['_']
    drop_cols_containing = ['@']
    # Conversion settings
    zero_repr = np.nan
    index_name = 'IS_genomicID'
    cols_for_index = ['chr', 'integration_locus', 'strand']
    # Func
    def parse_matrix (path):
        verbosePrint("> loading {path} ... ".format(path=str(path)))
        df = pd.read_csv(path,index_col=index_col, encoding=input_matrixes_encoding,sep=input_matrixes_sep,parse_dates=False)
        df.reset_index(inplace=True)
        return df
    def clean_matrix (df):
        verbosePrint("> detect columns to remove ...")
        drop_list = []
        for label in df.columns:
            if label in cols_to_drop:
                drop_list.append(label)
            else:
                for chars in drop_cols_starting_with:
                    if label.startswith(chars):
                        drop_list.append(label)
                for substring in drop_cols_containing:
                    if (substring in label):
                        drop_list.append(label)
        drop_list = list(set(drop_list))
        #verbosePrint("  drop_list: {drop_list}".format(drop_list=str(drop_list)))
        verbosePrint("> dropping columns ...")
        return df.drop(drop_list, axis=1)
    def convert (df):
        verbosePrint("> converting ...")
        new_index = df[cols_for_index].astype(str).apply(lambda x: 'chr'+'_'.join(x), axis=1)
        new_index.name = index_name
        return df.drop(cols_for_index, axis=1).set_index(new_index).replace(0, zero_repr)
    # Call
    verbosePrint("\n[LOAD MATRIX(ES) - Old Format]")
    matrixes = [convert(clean_matrix(parse_matrix(p))) for p in paths]
    verbosePrint("[DONE]")
    return matrixes

def export_matrixes_as_new (matrixes, paths):
    '''
    Purpose: export a Pandas DataFrame as text file
             NEW FORMAT
    
    IN: 
    matrixes - list of Pandas DataFrame to export
    path - list of paths exploited in writing
    
    OUT:
    0
    '''
    # Write settings
    output_matrixes_encoding = 'utf-8'
    output_matrix_sep = '\t'
    output_na_rep = ''
    # Func
    def write_df_as_new (df, path):
        verbosePrint("> writing {path} ... ".format(path=str(path)))
        df.to_csv(path_or_buf=path,
                  encoding=output_matrixes_encoding,
                  sep=output_matrix_sep,
                  na_rep=output_na_rep)
        verbosePrint("  done!")
    verbosePrint("\n[EXPORT MATRIX(ES) - As New Format]")
    [write_df_as_new(df, path) for df, path in zip (matrixes, paths)]
    verbosePrint("[DONE]")
    return 0
    
def export_matrixes_as_old (matrixes, paths):
    '''
    Purpose: export a Pandas DataFrame as text file
             OLD FORMAT
    IN: 
    matrixes - list of Pandas DataFrame to export
    path - list of paths exploited in writing
    
    OUT:
    0
    '''
    # Write settings
    output_matrixes_encoding = 'utf-8'
    output_matrix_sep = '\t'
    output_na_rep = ''
    # Conversion settings
    chr_col_name = 'chr'
    locus_col_name = 'integration_locus'
    strand_col_name = 'strand'
    total_col_name = 'all'
    # Func
    def convert (df):
        index = pd.Series(df.index)
        chr_col = index.str.split('_').str[0].apply(lambda x: x.replace('chr', ''))
        chr_col.name = chr_col_name
        locus_col = index.str.split('_').str[1]
        locus_col.name = locus_col_name
        strand_col = index.str.split('_').str[2]
        strand_col.name = strand_col_name
        all_col = pd.Series(df.sum(axis=1)).reset_index(drop=True)
        all_col.name = total_col_name
        return pd.concat([chr_col, locus_col, strand_col, df.reset_index(drop=True), all_col], axis=1)
    def write_df_as_old (df, path):
        verbosePrint("> processing {path} ...".format(path=str(path)))
        verbosePrint("  > converting ...")
        df_old = convert(df)
        verbosePrint("  > writing ...")
        df_old.to_csv(path_or_buf=path,
                  encoding=output_matrixes_encoding,
                  sep=output_matrix_sep,
                  na_rep=output_na_rep,
                  index=False)
        verbosePrint("    done!")
    verbosePrint("\n[EXPORT MATRIX(ES) - As Old Format]")
    [write_df_as_old(df, path) for df, path in zip (matrixes, paths)]
    verbosePrint("[DONE]")
    return 0

def check_and_build_outdir (outdir):
    '''
    '''
    def buildOUTDIR(ground_dir, *subfolders):
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
        # Check ground_dir: try to create it if not exists
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
        # Check ground_dir: write permissions (needed in case you don't give any *subfolders and ground_dir already exist)
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
            # Check: try to create OUTDIR if not exists
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
        # Check: write permissions (needed in case you specify *subfolders that already exists)
        if not os.access(OUTDIR, os.W_OK):
            print "\n[ERROR] You don't have write permissions in OUTDIR='{OUTDIR}'".format(OUTDIR=str(OUTDIR))
            sys.exit("\n[QUIT]\n")
        # Return OUTDIR
        verbosePrint("> OUTDIR: {OUTDIR}".format(OUTDIR=str(OUTDIR)))
        verbosePrint("> write permission: OK.")
        return OUTDIR
    
    if outdir is not None:
        items = outdir.split(os.sep)
        buildOUTDIR(items[0]+os.sep, *tuple(items[1:]))

def build_output_paths (outdir, suffix, *input_paths):
    '''
    Purpose: build a list of output_paths from *input_paths
    IN: 
    outdir - None or an absolute dir path. If None, output_paths will
             point to the same directory of the related input_path.
             Else, output_paths will point ALL to the same outdir
    suffix - a label + extension to avoid input overwrite
    
    OUT:
    list of output_paths
    '''
    # suffix embeds file extension
    def split_inpath (inpath):
        dirpath, name_and_ext = os.path.split(inpath)
        name, ext = os.path.splitext(name_and_ext)
        return dirpath, name, ext  # ext has dot
    def build_outpath (dirpath, name, suffix):
        return os.path.normpath(os.path.join(dirpath, name+suffix))
    verbosePrint("\n[CREATE AND CHECK OUT PATH(S)]")
    # list to return
    output_paths = []
    # user-specified outdir
    if outdir is not None:
        outdir = str(outdir)
        try:
            outdir = os.path.normpath(outdir)
        except Exception, err_message:
            print "\n[ERROR] outdir='{out_path}' is not formatted as valid dir path!".format(outdir=str(outdir))
            print "os.path.normpath returned: ", err_message
            sys.exit("\n[QUIT]\n")
        outdir = os.path.expanduser(outdir)
        outdir = os.path.expandvars(outdir)
        if not os.path.isabs(outdir):
            print "\n[ERROR] outdir must be an absolute dir path! Your outdir='{outdir}'".format(outdir=str(outdir))
            sys.exit("\n[QUIT]\n")
        check_and_build_outdir (outdir)
    # loop over input_paths to fill output_paths list
    for inpath in input_paths:
        dirpath, name, ext = split_inpath(inpath)
        if outdir is None:
            output_paths.append(build_outpath(dirpath, name, suffix))
        else:
            output_paths.append(build_outpath(outdir, name, suffix))
    verbosePrint("...OK!")
    verbosePrint("[DONE]")
    return output_paths  # list


### MAIN FUNCS ############################################################################################################################################################################

def main_old_to_new (*input_paths):
    verbosePrint("\n[START]")
    # manage entire directories along with single file paths  - returned input_paths is a list
    input_paths = input_paths_from_dir (*input_paths)
    # check and norm input paths - returned input_paths is a list
    check_input_paths (*input_paths)
    input_paths = normalize_input_paths (*input_paths)
    # build output_paths
    output_paths = build_output_paths (outdir, old_to_new_suffix, *input_paths)
    # load matrixes - returned matrixes is always a list of Pandas DataFrame(s), from now on
    matrixes = import_old_matrixes (*input_paths)
    # export matrix in new format
    export_matrixes_as_new (matrixes, output_paths)
    verbosePrint("\n[END]\n")
    return output_paths


def main_new_to_old (*input_paths):
    verbosePrint("\n[START]")
    # manage entire directories along with single file paths  - returned input_paths is a list
    input_paths = input_paths_from_dir (*input_paths)
    # check and norm input paths - returned input_paths is a list
    check_input_paths (*input_paths)
    input_paths = normalize_input_paths (*input_paths)
    # build output_paths
    output_paths = build_output_paths (outdir, new_to_old_suffix, *input_paths)
    # load matrixes - returned matrixes is always a list of Pandas DataFrame(s), from now on
    matrixes = import_new_matrixes (*input_paths)
    # export matrix in old format
    export_matrixes_as_old (matrixes, output_paths)
    verbosePrint("\n[END]\n")
    return output_paths

    
### MAIN #################################################################################################################################################################################

if __name__ == '__main__':
    
        if args.command == 'old_to_new':
            main_old_to_new (*input_paths)
        elif args.command == 'new_to_old':
            main_new_to_old (*input_paths)
        else:
            print "[CODING ERROR]"
            sys.exit("\n[QUIT]\n")


