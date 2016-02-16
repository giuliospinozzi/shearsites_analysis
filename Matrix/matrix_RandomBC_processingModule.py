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

#+++++++++++++++++++++++++++++++++++++++ FUNCTIONS ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

##############################################################################################################################################################################
### BUILD EXHAUSITEVE DATAFRAME FROM POOL_alldata_dict = matrix_randomBC_dataModule.loadDataFiles(...)                                                                       #
### NOTE: buildExhaustiveDataFrame will be improved to manage ALL data at the end. Up to now, is just like df = buildDataFrame() with 'header_list' column added             #
##############################################################################################################################################################################

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
    columns = ['barcode', 'genomic_coordinates', 'header', 'shearsite', 'r1_chr', 'r1_end', 'r1_quality', 'r1_start', 'r1_strand', 'r2_chr', 'r2_end', 'r2_quality', 'r2_start', 'r2_strand', 'randomBC']
    exhaustive_df = pd.DataFrame(l, columns=columns)
    # DONE! Here the structure is actually exhaustive
    
    # CREATE exhaustive_df JUST LIKE df = buildDataFrame() WITH 'header' COLUMN ADDED (HEADER LIST)
    # here:
    # 1 - drop columns
    # 2 - group-by (columns - ['header'] - columns_to_drop) = ['barcode', 'genomic_coordinates', 'shearsite', 'randomBC']
    # 3 - convert the only 'grouped column' (namely 'header') to list --> it's the header list now!
    # 4 - rename 'header' column as 'header_list'
    # 5 - reset index
    # 6 - add 'seq_count' column evaluating the len() of each list of headers
    columns_to_drop = ['r1_chr', 'r1_end', 'r1_quality', 'r1_start', 'r1_strand', 'r2_chr', 'r2_end', 'r2_quality', 'r2_start', 'r2_strand']
    exhaustive_df = exhaustive_df.drop(columns_to_drop, 1)
    columns_to_exclude_in_groupby = ['header'] + columns_to_drop
    columns_to_groupby = []
    for x in columns:
        if x not in columns_to_exclude_in_groupby:
            columns_to_groupby.append(x)
    exhaustive_df = exhaustive_df.groupby(columns_to_groupby)['header'].apply(lambda x: x.tolist())
    exhaustive_df = pd.DataFrame(exhaustive_df)
    exhaustive_df.rename(columns={'header': 'header_list'}, inplace=True)
    exhaustive_df = exhaustive_df.reset_index()
    # add seq_count column
    exhaustive_df['seq_count'] = exhaustive_df['header_list'].apply(lambda x: len(x))
    # return
    return exhaustive_df   # up to now: like df = buildDataFrame() with 'header_list' column added


##############################################################################################################################################################################
### BUILD DATAFRAME FROM POOL_IS_dict = matrix_randomBC_dataModule.loadDataFiles(...)                                                                                        #
### returned df is designed to embed the smallest data content to perform (almost) all the tasks available in this framework                                                 #
### NOTE: columns=['barcode', 'genomic_coordinates', 'shearsite', 'randomBC', 'seq_count'] is MANDATORY, ORDERED AND HIERARCHICAL                                            #
##############################################################################################################################################################################

def buildDataFrame(POOL_IS_dict):
    l=[]
    [l.append(list(k) + [v]) for k,v in flattenDict(POOL_IS_dict).items()]
    df = pd.DataFrame(l, columns=['barcode', 'genomic_coordinates', 'shearsite', 'randomBC', 'seq_count'])
    return df


##############################################################################################################################################################################
### BUILD MATRIX-LIKE DATAFRAMES (ISs as rows, SAMPLES as columns) FROM ANY_DF                                                                                               #
### NOTE: in the following, 'any_df' stands for df=buildDataFrame(...) as well as exhaustive_df=buildExhaustiveDataFrame(...)                                                #
##############################################################################################################################################################################

def buildSeqCountMatrix(any_df, sample_column='barcode'):
    # Check for required columns
    required_columns = set([sample_column, 'genomic_coordinates', 'seq_count'])
    if required_columns.intersection(set(any_df.columns)) != required_columns:
        print "\n[ERROR] buildSeqCountMatrix wrong input! Columns {required_columns} are required. Given dataframe has {columns_found}.".format(required_columns=str(humanSorted(list(required_columns))), columns_found=str(humanSorted(list(any_df))))
        sys.exit("\n[QUIT]\n")
    # pivot to get SEQUENCE COUNTS of ISs in df
    seqCount_matrix = pd.pivot_table(any_df, index='genomic_coordinates', columns=sample_column, values='seq_count', aggfunc='sum', margins=False)
    return seqCount_matrix

def buildShsCountMatrix(any_df, sample_column='barcode'):
    # Check for required columns
    required_columns = set([sample_column, 'genomic_coordinates', 'shearsite'])
    if required_columns.intersection(set(any_df.columns)) != required_columns:
        print "\n[ERROR] buildShsCountMatrix wrong input! Columns {required_columns} are required. Given dataframe has {columns_found}.".format(required_columns=str(humanSorted(list(required_columns))), columns_found=str(humanSorted(list(any_df))))
        sys.exit("\n[QUIT]\n")
    # pivot to get counts of distinct SHEARSITES for each ISs in df (unique method exploited by lambda in aggfunc is a method of pd.Series class)
    # NOTE: here MARGINS are not as expected (not partial 'sums' of course, but the 'overall' len(x.unique()) column!!)
    ShsCount_matrix = pd.pivot_table(any_df, index='genomic_coordinates', columns=sample_column, values='shearsite', aggfunc=lambda x: len(x.unique()), margins=False)
    return ShsCount_matrix
    
def buildBarcodeCountMatrix(any_df, sample_column='barcode'):
    # try to take required columns -> new DF
    required_columns = [sample_column, 'genomic_coordinates', 'randomBC']
    l = []
    try:
        [l.append(any_df.loc[:,c].to_frame()) for c in required_columns]
    except:
        print "\n[ERROR] buildBarcodeCountMatrix wrong input! Columns {required_columns} are required. Given dataframe has {columns_found}.".format(required_columns=str(humanSorted(required_columns)), columns_found=str(humanSorted(list(any_df))))
        sys.exit("\n[QUIT]\n")
    DF = pd.concat(l, axis=1, join='inner')
    # group-by remaining fields but randomBC - > duplicate entries may arise
    grouped = DF.groupby([sample_column,'genomic_coordinates'], as_index=False)
    # aggregate duplicates taking the 'count-distinct' as values
    barcodeCount_df = grouped.aggregate(lambda x: len(np.unique(x)))
    barcodeCount_df.rename(columns={'randomBC': 'randomBC_countDistinct'}, inplace=True)
    # pivot to shape barcodeCount_df as matrix of distinct RANDOM-BARCODES for each ISs in df
    barcodeCount_matrix = pd.pivot_table(barcodeCount_df, index='genomic_coordinates', columns=sample_column, values='randomBC_countDistinct')
    return barcodeCount_matrix
  
def buildCellCountMatrix(any_df, sample_column='barcode'):
    # Check for required columns
    required_columns = set([sample_column, 'genomic_coordinates', 'shearsite', 'randomBC'])  # 'shearsite' is redundant, just to be sure about the completeness of the hierarchy (exploited in pivot_table)
    if required_columns.intersection(set(any_df.columns)) != required_columns:
        print "\n[ERROR] buildShsCountMatrix wrong input! Columns {required_columns} are required. Given dataframe has {columns_found}.".format(required_columns=str(humanSorted(list(required_columns))), columns_found=str(humanSorted(list(any_df))))
        sys.exit("\n[QUIT]\n")
    # pivot to get counts of distinct SHEARSITE-RANDOMTAG couples for each ISs in df
    cellCount_matrix = pd.pivot_table(any_df, index='genomic_coordinates', columns=sample_column, values='randomBC', aggfunc=lambda x: len(x), margins=False)
    return cellCount_matrix

def buildFragmentEstimateMatrix(any_df, sample_column='barcode'):
    # import rpy2 and sonicLength package
    try:
        import rpy2.robjects as robjects
    except Exception, err_message:
        print "\n[ERROR] buildFragmentEstimateMatrix can't run due to some troubles with 'import rpy2.robjects as robjects'."
        print "'import rpy2.robjects as robjects' returned: ", err_message
        sys.exit("\n[QUIT]\n")
    try:
        from rpy2.robjects.packages import importr
    except Exception, err_message:
        print "\n[ERROR] buildFragmentEstimateMatrix can't run due to some troubles with 'from rpy2.robjects.packages import importr'."
        print "'from rpy2.robjects.packages import importr' returned: ", err_message
        sys.exit("\n[QUIT]\n")
    try:
        sonicLength = importr("sonicLength")
    except Exception, err_message:
        print "\n[ERROR] buildFragmentEstimateMatrix can't run due to some troubles occurred while loading 'sonicLength' R-library."
        print '''importr("sonicLength") returned: ''', err_message
        sys.exit("\n[QUIT]\n")
    # try to take required columns -> new DF
    required_columns = [sample_column, 'genomic_coordinates', 'shearsite']
    l = []
    try:
        [l.append(any_df.loc[:,c].to_frame()) for c in required_columns]
    except:
        print "\n[ERROR] buildFragmentEstimateMatrix wrong input! Columns {required_columns} are required. Given dataframe has {columns_found}.".format(required_columns=str(humanSorted(required_columns)), columns_found=str(humanSorted(list(any_df))))
        sys.exit("\n[QUIT]\n")
    DF = pd.concat(l, axis=1, join='inner')
    # detect samples
    sample_list = humanSorted(list(DF[sample_column].unique()))
    # list of dataframes with results
    results_list = []
    # loop over samples
    for sample in sample_list:
        # select sample data
        sample_DF = DF[DF[sample_column]==sample]
        sample_DF = sample_DF.drop(sample_column, 1)
        sample_DF = sample_DF.drop_duplicates()  #unique!
        locations_list, length_list = list(sample_DF.genomic_coordinates), list(sample_DF.shearsite)
        # Alias for sonicLength - estAbund calling
        estAbund = sonicLength.estAbund
        # Call estAbund and store returned object in results
        results = estAbund(robjects.StrVector(locations_list), robjects.FloatVector(length_list))
        theta = results.rx2("theta")
        estimations_theta = list(theta)
        locations_theta = list(theta.names)
        # shape data as dataframe: col=locations_theta, sample_column, estimations_theta
        tmp_dict = {'genomic_coordinates': locations_theta, sample_column:[sample]*len(estimations_theta), 'fragmentEstimate_count':estimations_theta}
        tmp_DF = pd.DataFrame(data=tmp_dict)
        # append results
        results_list.append(tmp_DF)
    # concat dataframes in results_list
    results_df = pd.concat(results_list)
    # pivot to shape results_df as matrix of ESTIMATES OF N PARENTAL FRAGMENT YIELDED BY SONICLENGTH for each ISs in df
    fragmentEstimate_matrix = pd.pivot_table(results_df, index='genomic_coordinates', columns=sample_column, values='fragmentEstimate_count')
    return fragmentEstimate_matrix


#++++++++++++++++++++++++++++++++++++++ MAIN and TEST +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

if __name__ == "__main__":
    pass