# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 11:29:56 2016

@author: stefano
"""

#++++++++++++++ Requested Package(s) Import +++++++++++++++#
import pandas as pd
import numpy as np
import collections
import matrix_RandomBC_globModule

#++++++++++++++++++++++ Global Vars +++++++++++++++++++++++#
# verbose = matrix_RandomBC_globModule.verbose


#++++++++++++++++++++++ Global Funcs ++++++++++++++++++++++#
# verbosePrint = matrix_RandomBC_globModule.verbosePrint
humanSorted = matrix_RandomBC_globModule.humanSorted

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
    seqCount_matrix = pd.pivot_table(df, index='genomic_coordinates', columns='barcode', values='seq_count', aggfunc='sum', margins=False)
    return seqCount_matrix

def buildShsCountMatrix(df):
    # pivot to get counts of distinct SHEARSITES for each ISs in df
    # unique method exploited by lambda in aggfunc is a method of pd.Series class
    ### NOTE: here MARGINS are not as expected (not partial 'sums' of course, but the 'overall' len(x.unique()) column!!)
    ShsCount_matrix = pd.pivot_table(df, index='genomic_coordinates', columns='barcode', values='shearsite', aggfunc=lambda x: len(x.unique()), margins=False)
    return ShsCount_matrix

def buildBarcodeCountMatrix(df):
    # drop useless data from df and get DF
    DF = df.drop('shearsite', 1)
    DF = DF.drop('seq_count', 1)
    # group-by remaining fields but randomBC - > duplicate entries may arise
    grouped = DF.groupby(['barcode','genomic_coordinates'], as_index=False)
    # aggregate duplicates taking the 'count-distinct' as values
    barcodeCount_df = grouped.aggregate(lambda x: len(np.unique(x)))
    barcodeCount_df.rename(columns={'randomBC': 'randomBC_countDistinct'}, inplace=True)
    # pivot to shape barcodeCount_df as matrix of distinct RANDOM-BARCODES for each ISs in df
    barcodeCount_matrix = pd.pivot_table(barcodeCount_df, index='genomic_coordinates', columns='barcode', values='randomBC_countDistinct')
    return barcodeCount_matrix

def buildCellCountMatrix(df):
    # pivot to get counts of distinct SHEARSITE-RANDOMTAG couples for each ISs in df
    cellCount_matrix = pd.pivot_table(df, index='genomic_coordinates', columns='barcode', values='seq_count', aggfunc=np.count_nonzero, margins=False)
    return cellCount_matrix

def buildFragmentEstimateMatrix(df):
    # import rpy2 and sonicLength package
    import rpy2.robjects as robjects
    from rpy2.robjects.packages import importr
    sonicLength = importr("sonicLength")
    # detect samples
    barcode_list = humanSorted(list(df['barcode'].unique()))
    # list of dataframes with results
    results_list = []
    # loop over samples
    for barcode in barcode_list:
        # select sample data
        DF = df[df['barcode']==barcode]
        DF = DF.drop('barcode', 1)
        DF = DF.drop('randomBC', 1)
        DF = DF.drop('seq_count', 1)
        DF = DF.drop_duplicates()  #unique!
        locations_list, length_list = list(DF.genomic_coordinates), list(DF.shearsite)
        # Alias for sonicLength - estAbund calling
        estAbund = sonicLength.estAbund
        # Call estAbund and store returned object in results
        results = estAbund(robjects.StrVector(locations_list), robjects.FloatVector(length_list))
        theta = results.rx2("theta")
        estimations_theta = list(theta)
        locations_theta = list(theta.names)
        # shape data as dataframe: col=locations_theta, barcode, estimations_theta
        tmp_dict = {'genomic_coordinates': locations_theta, 'barcode':[barcode]*len(estimations_theta), 'fragmentEstimate_count':estimations_theta}
        tmp_DF = pd.DataFrame(data=tmp_dict)
        # append results
        results_list.append(tmp_DF)
    # concat dataframes in results_list
    results_df = pd.concat(results_list)
    # pivot to shape results_df as matrix of ESTIMATES OF N PARENTAL FRAGMENT YIELDED BY SONICLENGTH for each ISs in df
    fragmentEstimate_matrix = pd.pivot_table(results_df, index='genomic_coordinates', columns='barcode', values='fragmentEstimate_count')
    return fragmentEstimate_matrix




#++++++++++++++++++++++++++++++++++++++ MAIN and TEST ++++++++++++++++++++++++++++++++++++++#

if __name__ == "__main__":
    pass