# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 11:29:56 2016

@author: stefano
"""

#++++++++++++++ Requested Package(s) Import +++++++++++++++#
import pandas as pd
import numpy as np
import collections

#++++++++++++++++++++++ Global Vars +++++++++++++++++++++++#
# import matrix_RandomBC_globModule
# verbose = matrix_RandomBC_globModule.verbose


#+++++++++++++++++++++++++++++++++++++++ FUNCTIONS +++++++++++++++++++++++++++++++++++++++++#

def flattenDict(d, parent_key='', sep='@@'):
    '''
    - sep is used only in computation, so a 'robust' separator is advised
    - parent_key can be used to add a fixed string on the key begins while flattening
    '''
    items = []
    for k, v in d.items():
        new_key = parent_key + sep + str(k) if parent_key else k
        if isinstance(v, collections.MutableMapping):
            items.extend(flattenDict(v, new_key, sep=sep).items())
        else:
            items.append((tuple(new_key.split(sep)), v))
    return dict(items)
    
def buildDataFrame(POOL_IS_dict):
    l=[]
    [l.append(list(k) + [v]) for k,v in flattenDict(POOL_IS_dict).items()]
    df = pd.DataFrame(l, columns=['barcode', 'genomic_coordinates', 'shearsite', 'randomBC', 'seq_count'])
    return df
    
def buildSeqCountMatrix(df):
    # pivot to get SEQUENCE COUNTS of ISs in df
    seqCount_matrix = pd.pivot_table(df, index='genomic_coordinates', columns='barcode', values='seq_count', aggfunc='sum', margins=True)
    return seqCount_matrix
    
def buildCellCountMatrix(df):
    # pivot to get counts of distinct SHEARSITE-RANDOMTAG couples for each ISs in df
    cellCount_matrix = pd.pivot_table(df, index='genomic_coordinates', columns='barcode', values='seq_count', aggfunc=np.count_nonzero, margins=True)
    return cellCount_matrix
    
def buildShsCountMatrix(df):
    # pivot to get counts of distinct SHEARSITES for each ISs in df
    # unique method exploited by lambda in aggfunc is a method of pd.Series class
    ### NOTE: here MARGINS are not as expected (not partial 'sums' of course, but the 'overall' len(x.unique()) column!!)
    ShsCount_matrix = pd.pivot_table(df, index='genomic_coordinates', columns='barcode', values='shearsite', aggfunc=lambda x: len(x.unique()))
    return ShsCount_matrix


#++++++++++++++++++++++++++++++++++++++ MAIN and TEST ++++++++++++++++++++++++++++++++++++++#

if __name__ == "__main__":
    pass