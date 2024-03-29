#!/usr/bin/python
# -*- coding: utf-8 -*-


"""
Created on Fri Aug  5 2016

@author: Stefano Brasca
"""


import sys
import os
import pandas as pd


__author__ = "Stefano Brasca"
__copyright__ = "SR-TIGET"
__credits__ = ["Stefano Brasca", "Andrea Calabria", "Giulio Spinozzi"]
__version__ = "1.0"
__maintainer__ = "Stefano Brasca"
__email__ = "brasca.stefano@hsr.it"
__status__ = "Testing"


### APPUNTI DI SVILUPPO ##############################################################################################################
### # COSE DA MIGLIORARE
### 1) tutti i 'sys.exit' e i 'print' associati sarebbero da sostituire con 'raise' e 'sys.stderr.write'
### # NUOVE FEATURES
### 1) Sarebbe bello introdurre un arg che specifica altre colonne da considerare come index (es annotazioni varie):
### il comportamento sarebbe consistenze perchè "spezzerebbe IS diversamente annotati" (possibili tricky usages).
### Una funzione per annotare matrici grezze completerebbe il quadro (uso come "riformattatore", a singolo *input_paths).
### Per rendere il programma più fruibile si possono pensare args speciali (tipo --input_format_CNR", "--human_annotation", ...)
### che se non dati fa niente, altrimenti fanno override di certi altri arg (se non uguali ai loro default, possibili tricky usages)
######################################################################################################################################


### NOTE ABOUT DEFAULT VALUES (as prog or module) #########
# These variables are hidden inside funcs (globally called)
# and they do not appear as expicit kwargs. The reason is
# that python 2.7 does not support the synthax:
#     func(*args, kwarg1='stuff', kwarg2='items', ...)
# while python 3 does.
# The only option could have been something like:
#     func(*args, kwargs**)
# and this could be useful for flexible usage, but 
# can be esily done anytime, by necessity, and does
# not improve neither readability nor maintainability
###########################################################


### NOTE ABOUT '\t' #######################################
# '\t' is very challanging to hanlde. Here, in exchange for
# code readability, all go smooth, nevertheless arg parse
# argument error shows '\\t' as an available option,
# instead of '\t'. This can't be fixed. However the help
# [-h] shows it correctly.
###########################################################


### DEFAULT VALUES AS MODULE ##############################################################################
# Print on screen
verbose = True
print_time = True  # eval if verbose is True
# About output
output_matrix_sep = '\t'
output_matrixes_encoding = 'utf-8'
output_na_rep = ''  # how to WRITE missing data
# About input
input_matrixes_sep = '\t'
input_matrixes_encoding = 'utf-8'
input_matrixes_prefix = None  # None or a list paired with input_paths

use_attributes = None  # None or list of int 0-based (do relabelling)
old_attibute_sep = '_'  # Exploited in relabelling (use_attributes is not None) but also by add_prefix!
new_attribute_sep = '_'  # It makes sense only if use_attributes is not None (do relabelling)
###########################################################################################################


if __name__ == '__main__':
    
    description = """
+---+ GENERAL DESCRIPTION +---------------------------------------------------+

This program takes as input one or more matrix files and unify them, row and
column wise, by summation. A single file is created as output.
With one input matrix only, the program acts similarly to a "re-formatter".

+---+ DETAILS +---------------------------------------------------------------+

INPATHs (required argument(s)) are one or more absolute paths of matrix files.

Matrix files are expected to be plain text files (encoding can be specified by
"--input_file_encoding" optional argument), tabulated with a file separator
("--input_file_separator" optional argument).

First column is expected to host row IDs (IS IDs or whatever) as well as first
row is interpreted as column labels. The cell in common is taken as
"row ID name" and can be any.

The rest of the matrix is taken as 'content' (no special 'margins' or whatever
are allowed). Matrix content of each cell must be numeric or empty.

INPUT COMPLIANCE IS UP TO THE USER AND NO CONTROLS ARE PERFORMED.

You can redefine column labels of input matrix, specifying which attributes you
wish to keep, with a sequence of indexes 0-based: you might want to keep just a
subset (columns-to-group-like behaviour), re-shuffle them or whatever, even
calling same attributes more than once ("--input_data_attributes" optional
argument"). Up to now, input/output attribute separators are statically set
to '_'.

In case of duplicated row IDs or column labels, the program merges them by
summation (mandatory in order to provide results consistent with expectations).

You can also add a signature to columns belonging to distinct matrix file,
before unifying data ("--input_data_prefix" optional argument). Obviously
prefixes impact on unification, that is column-label-based (you have to pay
attention in exchange for possible tricky usages).

OUTPATH (required argument) is the absolute path of the output matrix file.
Output file features (text encoding, field separator, empty cell content) can
be tuned through specific optional arguments ("--output_file_encoding",
"--output_file_separator", "--output_file_na_representation").

+---+ ABOUT ARGUMENTS: +------------------------------------------------------+ """
    
    ### DEFAULT VALUES AS PROGRAM ##################################################################################################
    # Print on screen
    verbose = False
    print_time = True  # eval if verbose is True
    # About output
    output_matrix_sep = '\t'  # PARSER HELP NEEDS TO BE MANUALLY UPDATED AFTER CHANGE
    output_matrixes_encoding = 'utf-8'
    output_na_rep = ''  # how to WRITE missing data   # PARSER HELP NEEDS TO BE MANUALLY UPDATED AFTER CHANGE
    # About input
    input_matrixes_sep = '\t'  # PARSER HELP NEEDS TO BE MANUALLY UPDATED AFTER CHANGE
    input_matrixes_encoding = 'utf-8'
    input_matrixes_prefix = None  # MUST BE NONE HERE. If parsed as argument, will be a list paired with input_paths
    use_attributes = None  # None or list/tuple of int 0-based, if parsed, for relabelling
    old_attibute_sep = '_'  # NOT PARSED UP TO NOW. Exploited in relabelling (use_attributes is not None) but also by add_prefix!
    new_attribute_sep = '_'  # NOT PARSED UP TO NOW. It makes sense only if use_attributes is not None (do relabelling)
    ################################################################################################################################
    
    ### ARGS ############################################################################################################################
    import argparse
    class MyParser(argparse.ArgumentParser): 
       def error(self, message):
          print 'error: %s\n' % message
          print ""
          self.print_help()
          sys.exit("\n[QUIT]\n")
    #parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)
    parser = MyParser(description=description, formatter_class=argparse.RawTextHelpFormatter)
    # positional
    parser.add_argument("outpath", metavar='OUTPATH', help="absolute file path for the output matrix")
    parser.add_argument("inpath", metavar='INPATHs', nargs='+', help="absolute file path(s) for input matrix(es)")
    # optional
    parser.add_argument("-v", "--verbose", action="store_true", help="verbose execution")
    # optional about output
    parser.add_argument("-out_sep", "--output_file_separator", metavar='FIELD_SEP', choices=[r'\t', ',', ';'], default=output_matrix_sep, help="field separator of the output matrix file.\nDefault is: {default}.\nAvailable option are [{tab}, ',', ';']".format(tab=r"'\t'", default=r"'\t'"))
    parser.add_argument("-out_enc", "--output_file_encoding", metavar='ENCODING', default=output_matrixes_encoding, help="encoding of the output matrix file.\nDefault is: '{default}'.\nAll the standard Python encodings are supported.\nE.g.: 'utf-8', 'ascii', ...".format(default=output_matrixes_encoding))
    parser.add_argument("-out_na", "--output_file_na_representation", metavar='NA', default=output_na_rep, help="how to write NA in the output matrix file.\nDefault is: {default}. E.g.: 0, '', ...".format(default="'' (empty)"))
    # optional about input
    parser.add_argument("-in_sep", "--input_file_separator", metavar='FIELD_SEP', choices=[r'\t', ',', ';'], default=input_matrixes_sep, help="field separator of input matrix file(s).\nDefault is: {default}.\nAvailable option are [{tab}, ',', ';']".format(tab=r"'\t'", default=r"'\t'"))
    parser.add_argument("-in_enc", "--input_file_encoding", metavar='ENCODING', default=input_matrixes_encoding, help="encoding of input matrix file(s).\nDefault is: '{default}'.\nAll the standard Python encodings are supported.\nE.g.: 'utf-8', 'ascii', ...".format(default=input_matrixes_encoding))
    parser.add_argument("-in_prefix", "--input_data_prefix", metavar='PREFIX', nargs='+', default=input_matrixes_prefix, help="allows to add a distinctive prefix to the columns\nbelonging to each input matrix.\nMust be a sequence of prefix(es) paired with\n'INPATHs' (empty prefixes, i.e. '', are allowed).".format(default=str(input_matrixes_prefix)))
    parser.add_argument("-in_attr", "--input_data_attributes", metavar='INDEXES', nargs='+', type=int, default=use_attributes, help="allows to take only some attributes from column\nlabels and re-compute the matrix(es). Must be a\nsequence of non-negative integers, indicating\nindexes of attributes to keep (0-based).\nDefault behaviour is 'keep all'.\nAttribute separator is assumed to be '{old_attibute_sep}' and\ncannot be set as arg in actual implementation.".format(old_attibute_sep=str(old_attibute_sep)))
    # Parse Args
    args = parser.parse_args()
    
    #####################################################################################################################################
    
    ### SET VALUES #########################################################################################################
    # positional
    out_path = args.outpath
    input_paths = args.inpath
    # optional
    if args.verbose != verbose:
        verbose = args.verbose
    # optional about output
    if str(args.output_file_separator).replace("\\t", "\t") != output_matrix_sep:
        output_matrix_sep = str(args.output_file_separator).replace("\\t", "\t")
    if args.output_file_encoding != output_matrixes_encoding:
        output_matrixes_encoding = args.output_file_encoding
    if args.output_file_na_representation != output_na_rep:
        output_na_rep = args.output_file_na_representation
    # optional about input
    if str(args.input_file_separator).replace("\\t", "\t") != input_matrixes_sep:
        input_matrixes_sep = str(args.input_file_separator).replace("\\t", "\t")
    if args.input_file_encoding != input_matrixes_encoding:
        input_matrixes_encoding = args.input_file_encoding
    if args.input_data_prefix != input_matrixes_prefix:
        if type(args.input_data_prefix) is list:
            if len(input_paths) != len(args.input_data_prefix):
                print "\n[ARGUMENT ERROR] -in_prefix/--input_data_prefix must be a sequence of prefix(es)\npaired with 'INPATHs'. See help.\n"
                #print "                  Details: len(input_paths)={i}!={p}=len(input_data_prefix).".format(i=str(len(input_paths)), p=str(len(args.input_data_prefix)))
                parser.print_help()
                sys.exit("\n[QUIT]\n")
        elif args.input_data_prefix is None:
            pass  #That's ok
        else:
            print "\n[CODING ERROR] something went wrong while handling -in_prefix/--input_data_prefix. If you used it explicitly, please retry without. If this is not your case, or if the problem persists, please report this bug to the author!"
            sys.exit("\n[QUIT]\n")
        input_matrixes_prefix = args.input_data_prefix  # list or None
    if args.input_data_attributes != use_attributes:
        use_attributes = args.input_data_attributes
    ########################################################################################################################


### FUNC #################################################################################################################################################################################

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

def pandas_is_updated(verbose=verbose):
    updated_version = '0.17.0'
    actual_version = str(pd.__version__)
    test = [updated_version, actual_version]
    if test[0] != humanSorted(test)[0]:
        s = "\n[Warning] This script was written to run under a more updated version of pandas package ({p} or later, your version is {n}).".format(p=str(updated_version), n=str(actual_version))
        verbosePrint(s, verbose=verbose)
        return False
    else:
        return True


def check_input_paths (*paths):
    '''
    Purpose: check input matrix file path(s);
             exit and explain why if something is wrong.
    
    IN: 
    *paths - complete path(s) of matrix file(s)
    
    OUT:
    0
    '''
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
        if not os.path.isfile(p):
            print "\n[ERROR] input path must point to an existing file! Your input path='{p}'".format(p=str(p))
            sys.exit("\n[QUIT]\n")
    verbosePrint("...OK!")
    verbosePrint("[DONE]")
    return 0

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

def import_matrixes (*paths):
    '''
    Purpose: load matrix file(s) as Pandas DataFrame(s)
    
    IN: 
    *paths - complete path(s) of matrix file(s)
    
    OUT:
    matrixes - list of DataFrame objects
    '''
    def parse_matrix (path):
        verbosePrint("> loading {path} ... ".format(path=str(path)))
        m = pd.read_csv(path, index_col=0, encoding=input_matrixes_encoding, sep=input_matrixes_sep)
        m = m.replace(0, pd.np.NaN)
        return m
    verbosePrint("\n[LOAD MATRIX(ES)]")
    matrixes = [parse_matrix(p) for p in paths]
    verbosePrint("[DONE]")
    return matrixes

def redefine_matrixes_col_labels (use_attributes, *matrixes):
    '''
    Purpose: 
    return *matrixes relabelled, exploiting a sub-set of attributes
    from original column labels, specified by 0-based indexes in use_attributes
    (duplicates may arise). Attributes can be just reordered, or taken twice,
    there are no constraints about this, just set use_attributes properly.
    
    IN:
    use_attributes - supposed to be None or a list of ints, further
    checks will be performed inplace if needed. As list of ints, it
    specifies indexes (0-based) of attributes in column labels
    to be exploited for relabelling. 
    *matrixes - ...
    
    OUT:
    matrixes - list of DataFrame objects
    '''
    def old_labels_as_tuple_list (df, old_attibute_sep=old_attibute_sep):
        return [tuple(label.split(old_attibute_sep)) for label in df.columns]
    
    def new_labels_as_tuple_list (df, use_attributes):
        return [ tuple([o[i] for i in use_attributes]) for o in old_labels_as_tuple_list (df) ]
    
    def build_relabelling_dict (df, use_attributes, new_attribute_sep=new_attribute_sep):
        '''
        Purpose: build a relabelling_dict suitable for
                 df.rename(columns=relabelling_dict)
        '''
        old_labels = df.columns
        new_labels = [new_attribute_sep.join(t) for t in new_labels_as_tuple_list (df, use_attributes)]
        relabelling_dict = {}
        for o, n in zip(old_labels, new_labels):
            relabelling_dict[o] = n
        return relabelling_dict
    
    def redefine_df_col_labels (df, use_attributes, **kwargs):
        '''
        Purpose: return a new df relabelled according
        to use_attributes rules.
        Kwargs are needed only to verbosePrint
        operations: if provided, must be both
        "n" and "tot" (in the spirit of
        "processing n out of tot ...")
        '''
        if kwargs:
            verbosePrint("> processing {i} of {l} ... ".format(i=str(kwargs['n']), l=str(kwargs['tot'])))
        relabelling_dict = build_relabelling_dict (df, use_attributes)
        return df.rename(columns=relabelling_dict)
    
    def check_before_relabelling (df, use_attributes):
        '''
        Purpose: check df and return False if something
        is not ok.
        '''
        # get splitted column labels
        labels_as_tuple_list = old_labels_as_tuple_list (df)
        # create a list of len of splitted column labels
        tuple_len_list = [len(t) for t in labels_as_tuple_list]
        # non homogeneous label structure (also due to old_attibute_sep)
        if len(set(tuple_len_list)) > 1:
            verbosePrint("[ERROR] non homogeneous label structure within input matrix(es)!")
            return False
        # single field labels or wrong old_attibute_sep
        if tuple_len_list[0] == 1:
            verbosePrint("[ERROR] can't find old_attibute_separator in (some) column labels!")
            return False
        # check label fields
        for t in labels_as_tuple_list:
            for f in t:
                if (f=='') or False:  # here further tests in place of 'False'
                    verbosePrint("[ERROR] bad split: old_attibute_separator generated empty attribute(s)!")
                    return False
        # use_attributes must be seq of in >= 0
        if min(use_attributes) < 0:
            verbosePrint("[ERROR] use_attributes schema must be sequence of int >= 0 !")
            return False
        # use_attributes goes beyond index of last field of column labels
        if tuple_len_list[0] <= max(use_attributes):
            verbosePrint("[ERROR] use_attributes schema goes beyond last index of column label attributes!")
            return False
        return True
    # Do nothing if use_attributes is None
    if use_attributes is None:
        return list(matrixes)  # return input *matrixes unchanged, as list
    # Here task to do
    verbosePrint("\n[REDEFINE MATRIX(ES) COLUMN LABELS]")
    verbosePrint("* use_attributes schema: {use_attributes}".format(use_attributes=str(tuple(use_attributes))))
    verbosePrint("* old_attibute_separator: {old_attibute_sep}".format(old_attibute_sep=str(old_attibute_sep)))
    verbosePrint("* new_attibute_separator: {new_attribute_sep}".format(new_attribute_sep=str(new_attribute_sep)))
    # Check
    for m in matrixes:
        check_m = check_before_relabelling (m, use_attributes)
        if check_m is False:
            verbosePrint("[SKIP THIS TASK]")
            return list(matrixes)  # return input *matrixes unchanged, as list
    # Relabel if check went straight
    matrixes = [redefine_df_col_labels (m, use_attributes, n=i, tot=len(matrixes)) for i,m in enumerate(matrixes, start=1)]
    verbosePrint("[DONE]")
    return matrixes

def merge_rows_by_index (df):
    '''
    Purpose: sum-up rows with same index
    Tested in version '0.15.0', works as expected
    (NaN taken as 0 for sum() purposes, but sum of
    NaN yields NaN)
    '''
    return df.groupby(df.index).sum()
    
def merge_cols_by_label (df):
    '''
    Purpose: sum-up colss with same label
    Tested in version '0.15.0', works as expected
    (NaN taken as 0 for sum() purposes, but sum of
    NaN yields NaN)
    '''
    return df.groupby(df.columns, axis=1).sum()

def compact_df (df, **kwargs):
    '''
    Purpose: apply merge_rows_by_index and
    merge_cols_by_label only if needed by df.
    Kwargs are needed only to verbosePrint
    operations: if provided, must be both
    "n" and "tot" (in the spirit of
    "processing n out of tot ...")
    '''
    if kwargs:
        verbosePrint("> processing {i} of {l} ... ".format(i=str(kwargs['n']), l=str(kwargs['tot'])))
    if ((len(set(df.index)) != len(df.index)) or (len(set(df.columns)) != len(df.columns))):
        if len(set(df.index)) != len(df.index):
            if kwargs:
                verbosePrint("  > compacting rows ... ")
            df = merge_rows_by_index(df)
        if len(set(df.columns)) != len(df.columns):
            if kwargs:
                verbosePrint("  > compacting cols ... ")
            df = merge_cols_by_label(df)
    else:
        verbosePrint("  nothing to do. Skip!")
    return df

def compact_matrixes (*matrixes):
    '''
    Purpose: if *matrixes have duplicate column or row
    IDs, related columns or rows will be merged by
    summation.
    
    IN: *matrixes
    OUT: matrixes - list of ...
    '''
    verbosePrint("\n[COMPACT MATRIX(ES)]")
    matrixes = [compact_df(m, n=i, tot=len(matrixes)) for i,m in enumerate(matrixes, start=1)]
    verbosePrint("[DONE]")
    return matrixes

def add_prefixes_to_matrixes(prefix_sequence, matrix_sequence):
    '''
    Purpose: add prefixes to matrixes, pair-wise (one prefix in prefix_sequence
    applied to all columns of one matrix in matrix_sequence). If prefix_sequence
    is None, matrix_sequence is returned unchanged (as list).
    
    IN:
    prefix_sequence - a sequence of item(s) that support __str__
    matrix_sequence - a sequence of Pandas DataFrame(s)
    
    OUT:
    matrixes - list of input Pandas DataFrame(s) of length min(len(prefix_sequence), len(matrix_sequence))
    
    NOTE: if input sequences have the same length, OUT is exhaustive.
    '''
    # check arguments
    from collections import Sequence
    if prefix_sequence is not None:
        if not ((isinstance(prefix_sequence, Sequence) and not isinstance(prefix_sequence, basestring)) and (isinstance(matrix_sequence, Sequence) and not isinstance(matrix_sequence, basestring))):
            print "\n[ERROR] add_prefixes_to_matrixes takes only sequences as arguments!"
            sys.exit("\n[QUIT]\n")
    if prefix_sequence is None:
        return list(matrix_sequence)  # return input matrix_sequence unchanged, as list
    # def func
    def add_prefix (prefix, *matrixes, **kwargs):
        '''
        Purpose: add a 'prefix' (the same) to cols of *matrixes (Pandas df)
                 Kwargs are needed only to verbosePrint operations: if provided,
                 must be both "n" and "tot" (in the spirit of "processing n out of tot ...")
        Return: a matrix in case of one *matrixes argument, else a list.
        '''
        if kwargs:
            verbosePrint("> processing {i} of {l} (prefix='{p}') ... ".format(i=str(kwargs['n']), l=str(kwargs['tot']), p=str(prefix)))
        pref_m = [df.add_prefix(prefix+old_attibute_sep) for df in matrixes]
        if len(pref_m) == 1:
            return pref_m[0]
        else:
            return pref_m
    # run func
    verbosePrint("\n[ADD PREFIX TO DATA]")
    matrixes = [add_prefix (*pm, n=i, tot=len(matrix_sequence)) for i,pm in enumerate(zip(prefix_sequence, matrix_sequence), start=1)]
    verbosePrint("[DONE]")      
    return matrixes

#def unify_matrixes (*matrixes):
#    '''
#    Purpose: sum-up Pandas DataFrame(s), row and column wise
#    
#    IN: 
#    *matrixes - Pandas DataFrame(s)
#    
#    OUT:
#    final_matrix -  resulting Pandas DataFrame
#    '''
#    if len(matrixes) > 1:
#        verbosePrint("\n[UNIFY MATRIXES]")
#        final_matrix = pd.DataFrame()
#        i=0
#        for m in matrixes:
#            i += 1
#            verbosePrint("> processing {i} of {l} ... ".format(i=str(i), l=str(len(matrixes))))
#            final_matrix = final_matrix.add(m,
#                                            fill_value=0)
#        verbosePrint("[DONE]")
#        return final_matrix
#    else:
#        return matrixes[0]

def unify_matrixes (*matrixes):
    '''
    Purpose: sum-up Pandas DataFrame(s), row and column wise
    
    IN: 
    *matrixes - Pandas DataFrame(s)
    
    OUT:
    final_matrix -  resulting Pandas DataFrame
    '''
    if len(matrixes) > 1:
        verbosePrint("\n[UNIFY MATRIXES]")
        verbosePrint("> unifying ...")
        final_matrix = None
        if pandas_is_updated(verbose=False):
            final_matrix = pd.concat(matrixes, copy=False)
        else:
            final_matrix = pd.concat(matrixes)  # unexpected keyword 'copy'
        verbosePrint("> consolidating ...")
        final_matrix = compact_df(final_matrix)
        verbosePrint("[DONE]")
        return final_matrix
    else:
        return matrixes[0]

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
    
def build_outpath (out_path):
    '''
    Purpose:
    1) check if out_path is ok, normpath, expanduser, expandvars.
    2) exploit buildOUTDIR to create desired folder's tree till the last directory.
    
    IN: 
    out_path - complete abs path of the file you want to write
    
    OUT:
    out_path - complete abs path of the file you want to write
    (checked, normed, expanded and existing till the last directory)
    '''
    verbosePrint("\n[CHECK AND CREATE OUT PATH]")
    out_path = str(out_path)
    try:
        out_path = os.path.normpath(out_path)
    except Exception, err_message:
        print "\n[ERROR] out_path='{out_path}' is not formatted as valid path!".format(out_path=str(out_path))
        print "os.path.normpath returned: ", err_message
        sys.exit("\n[QUIT]\n")
    out_path = os.path.expanduser(out_path)
    out_path = os.path.expandvars(out_path)
    if not os.path.isabs(out_path):
        print "\n[ERROR] out_path must be an absolute path! Your out_path='{out_path}'".format(out_path=str(out_path))
        sys.exit("\n[QUIT]\n")
    items = out_path.split(os.sep)
    buildOUTDIR(items[0]+os.sep, *tuple(items[1:-1]))
    verbosePrint("[DONE]")
    return str(out_path)

def export_matrix (matrix, path):
    '''
    Purpose: export a Pandas DataFrame as text file
    
    IN: 
    matrix - Pandas DataFrame to export
    path - complete path exploited in writing
    
    OUT:
    0
    '''
    verbosePrint("\n[EXPORT MATRIX]")
    verbosePrint("> processing ...")
    matrix.to_csv(path_or_buf=path,
                  encoding=output_matrixes_encoding,
                  sep=output_matrix_sep,
                  na_rep=output_na_rep)
    verbosePrint(">>> file created: {path}".format(path=str(path)))
    verbosePrint("[DONE]")
    return 0



### MAIN FUNC ############################################################################################################################################################################

def main(out_path, *input_paths):
    verbosePrint("\n[START]")
    # check pandas version and warn user
    pandas_is_updated(verbose=True)
    # check and norm input paths - returned input_paths is a list
    check_input_paths (*input_paths)
    input_paths = normalize_input_paths (*input_paths)
    # check and norm out path and create folder's tree if does not exist - returned out_path is a str
    out_path = build_outpath (out_path)
    # load matrixes - returned matrixes is always a list of Pandas DataFrame(s), from now on
    matrixes = import_matrixes (*input_paths)
    # redefine col label strcture according to use_attributes (acts only if use_attributes is not None)
    matrixes = redefine_matrixes_col_labels (use_attributes, *matrixes)
    # merge (sum) rows/cols with same IDs/Labels
    matrixes = compact_matrixes (*matrixes)
    # add prefixes (acts only if input_matrixes_prefix is not None)
    matrixes = add_prefixes_to_matrixes(input_matrixes_prefix, matrixes)
    # merge matrixes: acts only if len(matrixes) > 1, else matrixes[0] is returned - final_matrix is a Pandas DataFrame
    final_matrix = unify_matrixes (*matrixes)
    # export
    export_matrix (final_matrix, out_path)
    verbosePrint("\n[END]\n")
    return out_path



### MAIN #################################################################################################################################################################################

if __name__ == '__main__':
    
    ### LAUNCH MAIN FUNC ##########
    main(out_path, *input_paths)  #
    ###############################
    
    
    