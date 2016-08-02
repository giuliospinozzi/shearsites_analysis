# -*- coding: utf-8 -*-
"""
Created on Fri Jul 29 14:47:18 2016

@author: stefano
"""

#++++++++++++++++++++++++++++ Requested Package(s) Import +++++++++++++++++++++++++++++#
import pandas as pd
import sys


#+++++++++++++++++++++++++++++++++++++++ FUNCTIONS +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#


def loadRefactored_asDataframe (path, compression=None):
    sep = '\t'
    encoding = 'utf-8'
    header = 0
    index_col = 'prod_header'
    usecols = ['prod_header', 'prod_chr', 'prod_locus', 'prod_end', 'prod_strand', 'ref_associationid']
    return pd.read_csv(path, sep=sep, compression=compression, encoding=encoding, header=header, index_col=index_col, usecols=usecols)


def loadFasta_asDataframe (path, compression=None):
    # supported compression: None or 'gzip'
    def loadGzFile_asRowList (path, eol='\n'):
        import gzip
        file_content = None
        with gzip.open(path, 'r') as f:
            file_content = f.read()
        file_as_list = file_content.split(eol)
        return file_as_list
        
    def loadFile_asRowList (path, eol='\n'):
        file_content = None
        with open(path, 'r') as f:
            file_content = f.read()
            file_as_list = file_content.split(eol)
        return file_as_list
        
    def getHeader (header_line):
        # header_line has '>' and a suffix after a space to be removed
        header = header_line[1:].split(' ')[0]
        return header
    
    def getSeq (seq_line):
        # here you can manipulate seq_line
        return seq_line
    
    # No path no Data!
    if path is None:
        return None
    # Get file_as_list (eol removed)
    file_as_list = None
    if compression is None:
        file_as_list = loadFile_asRowList(path)
    elif compression == 'gzip':
        file_as_list = loadGzFile_asRowList(path)
    else:
        print "\n[ERROR] unsupported compression={compression} for file={path}".format(compression=str(compression), path=str(path))
        sys.exit("\n[QUIT]\n")
    # Get tuple_list [(header0, seq0), (header1, seq1), ...]
    from itertools import izip
    tuple_list = [(getHeader(header_line), getSeq(seq_line)) for header_line, seq_line in izip(*[iter(file_as_list)]*2)]
    # Convert list to np array, create Dataframe and return
    import numpy as np
    dt = np.dtype([('header', np.str_, 100), ('randomBC', np.str_, 12)])
    tuple_list = np.array(tuple_list, dtype=dt)
    # returnd DF has seq headers as index and a column named 'randomBC'
    return pd.DataFrame.from_records(tuple_list, index='header')


def buildDataFrame (refactored_DF, rBC_fasta_DF=None, drop_headers=True):
    
    def join (refactored_DF, rBC_fasta_DF):
        raw_DF = None
        if rBC_fasta_DF is not None:
            raw_DF = refactored_DF.join(rBC_fasta_DF, how='left')
        else:
            raw_DF = refactored_DF.copy()
        return raw_DF
        # raw_DF is just a copy of refactored_DF with or without 'randomBC'
        
    def transform (raw_DF, drop_headers):
        ### set up new columns ###
        # create 'shearsite'
        raw_DF['shearsite'] = raw_DF[['prod_locus', 'prod_end']].apply(lambda x: abs(x[0]-x[1]), axis=1)
        # drop useless stuff
        raw_DF.drop('prod_end', axis=1, inplace=True)
        # create 'genomic_coordinates'
        raw_DF['genomic_coordinates'] = raw_DF[['prod_chr', 'prod_locus', 'prod_strand']].apply(lambda x: "_".join((x[0], str(x[1]), x[2])), axis=1)
        # drop useless stuff
        raw_DF.drop(['prod_chr', 'prod_locus', 'prod_strand'], axis=1, inplace=True)
        # rename ref_associationid as 'association_ID'
        raw_DF.rename(columns={'ref_associationid': 'association_ID'}, inplace=True)
        # reset index / keep headers
        raw_DF.reset_index(inplace=True, drop=drop_headers)
        # change columns order and return
        col = ['association_ID', 'genomic_coordinates', 'shearsite', 'randomBC', 'header']
        if 'randomBC' not in raw_DF.columns:
            col.remove('randomBC')
        if drop_headers:
            col.remove('header')
        long_DF = raw_DF[col]
        return long_DF
        # NOTE:
        #      0) long_DF has one row per read
        #      1) for sure long_DF has columns: ['association_ID', 'genomic_coordinates', 'shearsite']
        #      2) 'header' column may exist or not (arg drop_headers: True o False)
        #      3) 'randomBC' column may exist or not (dipende dall'input, v. join function)
        #      4) seq_count column does not exist yet (should be identically 1!)
        # long_DF is just a pointer to raw_DF: raw_DF is MODIFIED INPLACE along the function.
        
    def compact (long_DF):
        import numpy as np
        # set up grouping rule (all but 'header' if present)
        grouping = list(long_DF.columns)
        if 'header' in grouping:
            grouping.remove('header')
        # exhausitve df
        if 'header' in long_DF.columns:
            # create groups
            # work-around "the function does not reduce" error. This way I can get lists of headers in a new column named 'header_list'
            exhaustive_df = long_DF.groupby(grouping)['header'].apply(lambda x: x.tolist())
            exhaustive_df = pd.DataFrame(exhaustive_df)
            exhaustive_df.rename(columns={'header': 'header_list'}, inplace=True)
            exhaustive_df.reset_index(inplace=True)
            # add 'seq_count' column as len(header list)
            exhaustive_df['seq_count'] = exhaustive_df['header_list'].apply(lambda x: len(x))
            return exhaustive_df
        # standard df
        else:
            df = long_DF.copy()
            # add new seq_count column (identically 1 for each row)
            df['seq_count'] = [1]*len(df)
            # create groups
            grouped = df.groupby(grouping, as_index=False)
            # aggregate by summing up seq_count
            df = grouped.agg({'seq_count': np.sum})
            df.reset_index(inplace=True, drop=True)
            return df
    
    raw_DF = join(refactored_DF, rBC_fasta_DF)  # return a copy, inputs (refactored_DF and rBC_fasta_DF) are not modified
    long_DF = transform(raw_DF, drop_headers)   # raw_DF is modified inplace: long_DF is just a pointer to it
    any_df = compact(long_DF)  # long_DF is not modified, df is a new dataframe
    
    return any_df   # any_df can be of kind 'df' or 'exhaustive_df', with or without randomBC (kwargs drop_headers and rBC_fasta_DF)


#++++++++++++++++++++++++++++++++++++++ MAIN and TEST +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#


if __name__ == "__main__":
    
    ### Test vars ###
    refactored_path = '/home/stefano/Desktop/RandomBC_matrix_development/test_input/prova_refactored.tsv'
    fasta_path = '/home/stefano/Desktop/RandomBC_matrix_development/test_input/prova_fasta.fasta'
    ### Test code ###
    DF_refactored = loadRefactored_asDataframe(refactored_path)
    DF_randomBC = loadFasta_asDataframe(fasta_path)
    # do stuff for test purposes
    DF_refactored = DF_refactored.append(DF_refactored) #double data
    DF_refactored.index = range(1000, 1000+len(DF_refactored)) #forced tricky indexing to let something match
    DF_refactored.index.name = 'header' #fix index name
    DF_randomBC = DF_randomBC.append(DF_randomBC) #double data
    DF_randomBC.index = range(1000, 1000+len(DF_randomBC)) #forced tricky indexing to let something match
    DF_randomBC.index.name = 'header' #fix index name
    # Test buildDataFrame
    exhaustive_df = buildDataFrame (DF_refactored, rBC_fasta_DF=DF_randomBC, drop_headers=False)
    exhaustive_df_noTag = buildDataFrame (DF_refactored, rBC_fasta_DF=None, drop_headers=False)
    df = buildDataFrame (DF_refactored, rBC_fasta_DF=DF_randomBC, drop_headers=True)
    df_noTag = buildDataFrame (DF_refactored, rBC_fasta_DF=None, drop_headers=True)


