# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 11:29:56 2016

@author: stefano
"""

#++++++++++++++ Requested Package(s) Import +++++++++++++++#
import sys
import pandas as pd
import numpy as np
import matrix_RandomBC_globModule

#++++++++++++++++++++++ Global Vars +++++++++++++++++++++++#
# verbose = matrix_RandomBC_globModule.verbose


#++++++++++++++++++++++ Global Funcs ++++++++++++++++++++++#
# verbosePrint = matrix_RandomBC_globModule.verbosePrint
humanSorted = matrix_RandomBC_globModule.humanSorted
flattenDict = matrix_RandomBC_globModule.flattenDict

#+++++++++++++++++++++++++++++++++++++++ FUNCTIONS +++++++++++++++++++++++++++++++++++++++++#

  
def buildExhaustiveDataFrame(POOL_alldata_dict):
    
    # Prepare POOL_alldata_dict for flattening -> d
    d = {}
    for barcode in POOL_alldata_dict.keys():
        if barcode not in d.keys():
            d[barcode] = {}
        for IS_id in POOL_alldata_dict[barcode].keys():
            if IS_id not in d[barcode].keys():
                d[barcode][IS_id] = {}
            for header in POOL_alldata_dict[barcode][IS_id].keys():
                if header not in d[barcode][IS_id].keys():
                    d[barcode][IS_id][header] = [str(POOL_alldata_dict[barcode][IS_id][header]['length']),
                                                 POOL_alldata_dict[barcode][IS_id][header]['r1_chr'],
                                                 str(POOL_alldata_dict[barcode][IS_id][header]['r1_end']),
                                                 str(POOL_alldata_dict[barcode][IS_id][header]['r1_quality']),
                                                 str(POOL_alldata_dict[barcode][IS_id][header]['r1_start']),
                                                 POOL_alldata_dict[barcode][IS_id][header]['r1_strand'],
                                                 POOL_alldata_dict[barcode][IS_id][header]['r2_chr'],
                                                 str(POOL_alldata_dict[barcode][IS_id][header]['r2_end']),
                                                 str(POOL_alldata_dict[barcode][IS_id][header]['r2_quality']),
                                                 str(POOL_alldata_dict[barcode][IS_id][header]['r2_start']),
                                                 POOL_alldata_dict[barcode][IS_id][header]['r2_strand'],
                                                 POOL_alldata_dict[barcode][IS_id][header]['randomBC'],
                                                ]
                    # d[barcode][IS_id][header] = ['length', 'r1_chr', 'r1_end', 'r1_quality', 'r1_start', 'r1_strand', 'r2_chr', 'r2_end', 'r2_quality', 'r2_start', 'r2_strand', 'randomBC']
                else:
                    print "\n[ERROR] Duplicate header found in {barcode} - {IS_id} : '{header}'".format(barcode=barcode, IS_id=IS_id, header=header)
                    sys.exit("\n[QUIT]\n")
    
    # flatten d and build exhaustive_df
    l=[]
    [l.append(list(k) + v) for k,v in flattenDict(d).items()]
    columns = ['barcode', 'genomic_coordinates', 'header', 'length', 'r1_chr', 'r1_end', 'r1_quality', 'r1_start', 'r1_strand', 'r2_chr', 'r2_end', 'r2_quality', 'r2_start', 'r2_strand', 'randomBC']
    exhaustive_df = pd.DataFrame(l, columns=columns)
    
    # DONE! Here the structure is actually exhaustive
    
    # here drop columns and append headers through group-by
    exhaustive_df = exhaustive_df.drop(['r1_chr', 'r1_end', 'r1_quality', 'r1_start', 'r1_strand', 'r2_chr', 'r2_end', 'r2_quality', 'r2_start', 'r2_strand'], 1)
    for x in ['header', 'r1_chr', 'r1_end', 'r1_quality', 'r1_start', 'r1_strand', 'r2_chr', 'r2_end', 'r2_quality', 'r2_start', 'r2_strand']:
        columns.remove(x)
    exhaustive_df = exhaustive_df.groupby(columns)['header'].apply(lambda x: x.tolist())
    exhaustive_df = pd.DataFrame(exhaustive_df)
    exhaustive_df = exhaustive_df.reset_index()
    # add seq_count column
    exhaustive_df['seq_count'] = exhaustive_df['header'].apply(lambda x: len(x))
    
    return exhaustive_df   # up to now: like df = buildDataFrame() with 'header' column added (header list)

##############################################################################################################################################################################
### NOTE: buildExhaustiveDataFrame will be improved to manage ALL data at the end. Up to now, is just like df = buildDataFrame() with 'header' column added (header list)    #
###       Despite the structure is ongoing, buildXxxxMatrix below can be adapted right now!                                                                                  #
###       Functions exploiting pivot_table should be always ok (check it)                                                                                                    #
###       Conversely, function exploiting 'drop-and-groupby' should be converted in a 'extractcolumns-and-groupby' fashion!                                                  #
##############################################################################################################################################################################

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