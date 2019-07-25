#!/usr/bin/python3
# -*- coding: utf-8 -*-


"""
Created on Wed July 17 2019

@author: Adriano De Marino
"""


from __future__ import print_function
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
from functools import reduce
import os, sys, argparse, time, multiprocessing, math, itertools, datetime, re
from multiprocessing import Queue,Process
import warnings, gc
warnings.simplefilter(action='ignore', category=FutureWarning)


__author__ = ["Adriano De Marino","Stefano Brasca"]
__copyright__ = "SR-TIGET"
__credits__ = ["Adriano De Marino","Stefano Brasca", "Andrea Calabria", "Giulio Spinozzi"]
__version__ = "2.0"
__updates__ = "memory improving"
__maintainer__ = "Adriano De Marino"
__email__ = "demarino.adriano@hsr.it"
__status__ = "Testing"


### APPUNTI DI SVILUPPO ##############################################################################################################
### # COSE DA MIGLIORARE
### 1) tutti i 'sys.exit' e i 'print(' associati sarebbero da sostituire con 'raise' e 'sys.stderr.write'
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
g = '\033[92m'
e = '\033[0m'
w = '\033[93m'
r = '\033[91m'
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
          print( 'error: %s\n' % message)
          print( "")
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
                print( "\n[ARGUMENT ERROR] -in_prefix/--input_data_prefix must be a sequence of prefix(es)\npaired with 'INPATHs'. See help.\n")
                #print( "                  Details: len(input_paths)={i}!={p}=len(input_data_prefix).".format(i=str(len(input_paths)), p=str(len(args.input_data_prefix)))
                parser.print_help()
                sys.exit("\n[QUIT]\n")
        elif args.input_data_prefix is None:
            pass  #That's ok
        else:
            print( "\n[CODING ERROR] something went wrong while handling -in_prefix/--input_data_prefix. If you used it explicitly, please retry without. If this is not your case, or if the problem persists, please report this bug to the author!")
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
    Purpose: print( immediately if verbose is True, with
             date and time if print(_time is True too.
    
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
            nnl = y.count("\n", 0, int(len(x)/2+1))
            y = y.replace("\n", "", nnl)
            nl_str = "".join(['\n']*nnl)
            print(nl_str+"[{time}] ".format(time=strftime("%Y-%m-%d %H:%M:%S", localtime())), y)
        else:
            print( x)
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
        except( Exception, err_message):
            print( "\n[ERROR] input path='{p}' is not formatted as valid path!".format(p=str(p)))
            print( "os.path.normpath returned: ", err_message)
            sys.exit("\n[QUIT]\n")
        p = os.path.expanduser(p)
        p = os.path.expandvars(p)
        # check if path is absolute, read permission and if path points to an existing file
        if not os.path.isabs(p):
            print( "\n[ERROR] input path must be an absolute path! Your input path='{p}'".format(p=str(p)))
            sys.exit("\n[QUIT]\n")
        if not os.access(os.path.dirname(p), os.R_OK):
            print( "\n[ERROR] You have not read permission in folder='{f}'".format(f=str(os.path.dirname(p))))
            sys.exit("\n[QUIT]\n")
        if not os.path.isfile(p):
            print( "\n[ERROR] input path must point to an existing file! Your input path='{p}'".format(p=str(p)))
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

class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout


def is_categorical(array_like):
    return array_like.dtype.name == 'category'

def reduce_mem_usage(props):
    # start_mem_usg = props.memory_usage().sum() / 1024**2 
    # print("Memory usage of properties dataframe is: {:03.5f} MB".format(start_mem_usg))
    NAlist = [] # Keeps track of columns that have missing values filled in. 
    for col in props.columns:
        if is_categorical(props[col]):
            continue
        if np.issubdtype(props[col].dtype, np.datetime64):
            continue
        if props[col].dtype != object:  # Exclude strings
            
            # Print current column type
            # print("******************************")
            # print("Column: ",col)
            # print("dtype before: ",props[col].dtype)
            
            # make variables for Int, max and min
            IsInt = False
            mx = props[col].max()
            mn = props[col].min()
            
            # Integer does not support NA, therefore, NA needs to be filled
            if not np.isfinite(props[col]).all(): 
                NAlist.append(col)
                props[col].fillna(mn-1,inplace=True)  
                   
            # test if column can be converted to an integer
            asint = props[col].fillna(0).astype(np.int64)
            result = (props[col] - asint)
            result = result.sum()
            if result > -0.01 and result < 0.01:
                IsInt = True
            
            # Make Integer/unsigned Integer datatypes
            if IsInt:
                if mn >= 0:
                    if mx < 255:
                        props.loc[:,col] = props[col].astype(np.uint8)
                    elif mx < 65535:
                        props.loc[:,col] = props[col].astype(np.uint16)
                    elif mx < 4294967295:
                        props.loc[:,col] = props[col].astype(np.uint32)
                    else:
                        props.loc[:,col] = props[col].astype(np.uint64)
                else:
                    if mn > np.iinfo(np.int8).min and mx < np.iinfo(np.int8).max:
                        props.loc[:,col] = props[col].astype(np.int8)
                    elif mn > np.iinfo(np.int16).min and mx < np.iinfo(np.int16).max:
                        props.loc[:,col] = props[col].astype(np.int16)
                    elif mn > np.iinfo(np.int32).min and mx < np.iinfo(np.int32).max:
                        props.loc[:,col] = props[col].astype(np.int32)
                    elif mn > np.iinfo(np.int64).min and mx < np.iinfo(np.int64).max:
                        props.loc[:,col] = props[col].astype(np.int64)    
            
            # Make float datatypes 32 bit
            else:
                props.loc[:,col] = props[col].astype(np.float32)
            
            # Print new column type
            # print("dtype after: ",props[col].dtype)
            # print("******************************")
        else:
            num_unique_values = len(props[col].unique())
            num_total_values = len(props[col])
            if num_unique_values / num_total_values < 0.5:
                props.loc[:,col] = props[col].astype('category')
            else:
                continue
    # Print final result
    # print("___MEMORY USAGE AFTER COMPLETION:___")
    # mem_usg = props.memory_usage().sum() / 1024**2 
    # print("Memory usage is: {:03.5f} MB".format(mem_usg))
    # print("Saved {:03.2f} % memory usage".format((start_mem_usg-mem_usg)/start_mem_usg*100))
    return props, NAlist
    
def chunks_columns(l, n):
    return [l[l.columns[i:i+int(n)]] for i in range(0, len(l.columns), int(n))]

def reducer(job_id,data_slice,return_dict,tail):
    #are they share somethings?
    gb = data_slice.copy()
    # print(np.shape(gb))
    with HiddenPrints():
        return_dict[job_id] = reduce_mem_usage(gb)[0]
    tail.put(job_id)

def mapper(data, job_number,mem_usage_list = []):
    total = len(data.columns)
    chunk_size = total / job_number
    slice = chunks_columns(data, chunk_size)
    jobs = []
    manager = multiprocessing.Manager()
    return_dict = manager.dict()
    tail = manager.Queue()

    for i, s in enumerate(slice):
        j = multiprocessing.Process(target=reducer, args=(i, s, return_dict,tail))
        jobs.append(j)
        #mem_usage_list.append(process.memory_info().rss/1024**2)
    for j in jobs:
        j.start()
    for p in jobs:
        p.join()

    return return_dict.values()#,np.sum(mem_usage_list)


def import_matrixes (*paths):
    '''
    Purpose: load matrix file(s) as Pandas DataFrame(s)
    
    IN: 
    *paths - complete path(s) of matrix file(s)
    
    OUT:
    matrixes - list of DataFrame objects
    '''
    def parse_matrix(path):
        jobs=48
        verbosePrint("> loading {path} ... ".format(path=str(path)))
        chunksize = 100000
        df = pd.DataFrame()
        memory_total = []
        memory_final = []
        for CHUNK in pd.read_csv(path,sep=input_matrixes_sep,index_col=0,encoding=input_matrixes_encoding,chunksize=chunksize):
            memory_total.append(CHUNK.memory_usage().sum() / 1024**2)
            CHUNK = CHUNK.fillna(0)
            if len(CHUNK) > jobs:
                new = pd.concat(mapper(CHUNK,jobs),axis=1,sort=False)
            else:
                new = pd.concat(mapper(CHUNK,1),axis=1,sort=False)
            memory_final.append(new.memory_usage().sum() / 1024**2)
            df = df.append(new)
        verbosePrint("> Memory usage of properties dataframe is: {mb:03.3f} MB ".format(mb=np.sum(memory_total)))
        verbosePrint("> ___MEMORY USAGE AFTER COMPLETION:______: {ma:03.3f} MB".format(ma=np.sum(memory_final)))
        return df
        # m = pd.read_csv(path, index_col=0, encoding=input_matrixes_encoding, sep=input_matrixes_sep) old code
        # m = m.replace(0, pd.np.NaN) old code
        # return m old code
    verbosePrint("\n[LOAD MATRIX(ES)]")
    matrixes = [parse_matrix(p) for p in paths]
    verbosePrint("[DONE]")
    return matrixes

def updateInPlace(a,b):
    a.update(b)
    return a

def unify_n(s,super):
    
    piece = s.loc[:, (s != 0).any(axis=0)]
    piece_dict = piece.to_dict()
    for k,l in piece_dict.items():
        for cy in l:
            if l[cy] != 0:
                super.append((k+'..'+cy, l[cy]))

def merge_dkeys(data, job_number):
        
    total = len(data)
    chunk_size = total / job_number
    slice = chunks(data, chunk_size)
    jobs = []
    manager = multiprocessing.Manager()
    super = manager.list()

    for s in slice:
        j = multiprocessing.Process(target=unify_n, args=(s, super))
        jobs.append(j)
    for j in jobs:
        j.start()
    for p in jobs:
        p.join()

    return super

def unify_matrixes (*matrixes):
    '''
    Purpose: sum-up Pandas DataFrame(s), row and column wise
    
    IN: 
    *matrixes - Pandas DataFrame(s)
    
    OUT:
    final_matrix -  resulting Pandas DataFrame
    ADRIANO OPTIMISATION 
    '''
    out_path2 = out_path.split('.')[0]+'_TMP.'+out_path.split('.')[1]
    number_matrixes = len(matrixes)
    verbosePrint("\n[UNIFY MATRIXES]")
    verbosePrint("> unifying ...")
    final_matrixes = []

    for mat in matrixes:
        final_matrixes.append(merge_dkeys(mat, 56))

    c = reduce(updateInPlace, (Counter(dict(x)) for x in final_matrixes))
    sample = c.keys()
    value = c.values()

    samples = pd.Series([x.split('..')[0] for x in sample],name='nome').to_frame()
    valore = pd.Series([x for x in value],name='value').to_frame()
    IS = pd.Series([x.split('..')[1] for x in sample],name='IS').to_frame()

    matrix_final = samples.join([IS,valore])
    verbosePrint("[DONE]")

    return matrix_final

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
    except(Exception, err_message):
        print( "\n[ERROR] ground_dir='{ground_dir}' is not formatted as valid path!".format(ground_dir=str(ground_dir)))
        print( "os.path.normpath returned: ", err_message)
        sys.exit("\n[QUIT]\n")
    if os.path.isfile(OUTDIR):
        print( "\n[ERROR] ground_dir='{OUTDIR}' must be a folder path, not a file!".format(OUTDIR=str(OUTDIR)))
        sys.exit("\n[QUIT]\n")
    if not os.path.isabs(OUTDIR):
        print( "\n[ERROR] ground_dir='{OUTDIR}' must be an absolute path!".format(OUTDIR=str(OUTDIR)))
        sys.exit("\n[QUIT]\n")
    # Check ground_dir: try to create it if not exists
    if not os.path.exists(OUTDIR):
        try:
            os.makedirs(OUTDIR)
        except (Exception, err_message):
            print( "\n[ERROR] ground_dir='{OUTDIR}' is not a valid path or you don't have permissions to create it!".format(OUTDIR=str(OUTDIR)))
            print( "os.makedirs returned: ", err_message)
            sys.exit("\n[QUIT]\n")
        verbosePrint("> ground_dir created: {OUTDIR}".format(OUTDIR=str(OUTDIR)))
    else:
        verbosePrint("> ground_dir found: {OUTDIR}".format(OUTDIR=str(OUTDIR)))
    # Check ground_dir: write permissions (needed in case you don't give any *subfolders and ground_dir already exist)
    if not subfolders:
        if not os.access(OUTDIR, os.W_OK):
            print( "\n[ERROR] You don't have write permissions in ground_dir='{OUTDIR}'".format(OUTDIR=str(OUTDIR)))
            sys.exit("\n[QUIT]\n")
    # Loop over *subfolders
    for sf in subfolders:
        # Try OUTDIR=join(OUTDIR, sf) and formal check
        try:
            OUTDIR = os.path.normpath(os.path.join(OUTDIR, sf))
        except( Exception, err_message):
            print( "\n[ERROR] Cannot join OUTDIR='{OUTDIR}' and SUBFOLDER='{SUBFOLDER}' as a valid path!".format(OUTDIR=str(OUTDIR), SUBFOLDER=str(sf)))
            print( "os.path.normpath(os.path.join(...)) returned: ", err_message)
            sys.exit("\n[QUIT]\n")
        if os.path.isfile(OUTDIR):
            print( "\n[ERROR] *subfolders args must be folder(s), not file(s)! Input subfolders: {subfolders}.".format(subfolders=str(subfolders)))
            sys.exit("\n[QUIT]\n")
        # Check: try to create OUTDIR if not exists
        if not os.path.exists(OUTDIR):
            try:
                os.makedirs(OUTDIR)
            except (Exception, err_message):
                print( "\n[ERROR] OUTDIR='{OUTDIR}' + SUBFOLDER='{SUBFOLDER}' is not a valid path or you don't have permissions to create it!".format(OUTDIR=str(OUTDIR), SUBFOLDER=str(sf)))
                print( "os.makedirs returned: ", err_message)
                sys.exit("\n[QUIT]\n")
            verbosePrint("> subfolder created: {OUTDIR}".format(OUTDIR=str(OUTDIR)))
        else:
            verbosePrint("> subfolder found: {OUTDIR}".format(OUTDIR=str(OUTDIR)))
    # Check: write permissions (needed in case you specify *subfolders that already exists)
    if not os.access(OUTDIR, os.W_OK):
        print( "\n[ERROR] You don't have write permissions in OUTDIR='{OUTDIR}'".format(OUTDIR=str(OUTDIR)))
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
    except (Exception, err_message):
        print( "\n[ERROR] out_path='{out_path}' is not formatted as valid path!".format(out_path=str(out_path)))
        print( "os.path.normpath returned: ", err_message)
        sys.exit("\n[QUIT]\n")
    out_path = os.path.expanduser(out_path)
    out_path = os.path.expandvars(out_path)
    if not os.path.isabs(out_path):
        print( "\n[ERROR] out_path must be an absolute path! Your out_path='{out_path}'".format(out_path=str(out_path)))
        sys.exit("\n[QUIT]\n")
    items = out_path.split(os.sep)
    buildOUTDIR(items[0]+os.sep, *tuple(items[1:-1]))
    verbosePrint("[DONE]")
    return str(out_path)

def chunks(l, n):
    return [l[i:i+int(n)] for i in range(0, len(l), int(n))]

def export_matrix (matrix, path):
    '''
    Purpose: export a Pandas DataFrame as text file
    
    IN: 
    matrix - Pandas DataFrame to export
    path - complete path exploited in writing
    
    OUT:
    0
    '''
    # split a list into evenly sized chunks
    header = humanSorted(matrix.nome.unique().tolist())
    grp = [t for x,t in matrix.groupby('IS',sort=True)]

    def chunks(l, n):
        return [l[i:i+int(n)] for i in range(0, len(l), int(n))]

    verbosePrint(">>> file created: {path}".format(path=str(path)))
    
    pd.DataFrame(columns=header).to_csv(path_or_buf=path,
                  index_label='IS_genomicID',
                  encoding=output_matrixes_encoding,
                  sep=output_matrix_sep,
                  na_rep=output_na_rep)

    def exportazion(s,super,header):
        for gb in s:
            tmp = pd.concat([pd.DataFrame(columns=header),pd.pivot_table(gb,columns='nome',index='IS',values='value')],sort=False,copy=False)
            tmp.to_csv(path_or_buf=path, mode='a', header=False,sep='\t',na_rep=output_na_rep,encoding=output_matrixes_encoding)
    
    def mapper_expo(csv_list, job_number,header):
        copia = csv_list
        total = len(copia)
        chunk_size = total / job_number
        slice = chunks(copia, chunk_size)
        jobs = []
        manager = multiprocessing.Manager()
        super = manager.list()
        for s in slice:
            j = multiprocessing.Process(target=exportazion, args=(s, super, header))
            jobs.append(j)
        for j in jobs:
            j.start()
        for p in jobs:
            p.join()
        return super

    verbosePrint("\n[EXPORT MATRIX]")
    mapper_expo(grp,48,header) #slow, 48processor instead of 1 increas the speed of writing, but mix the order of the chr.
    verbosePrint("[DONE]")
    
    verbosePrint("\n[RE-LOADED MATRIX]")
    import_matrixes(path).sort_index(inplace=True).to_csv(path_or_buf=path,sep='\t',index_label='IS_genomicID',na_rep=output_na_rep,encoding=output_matrixes_encoding)
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

