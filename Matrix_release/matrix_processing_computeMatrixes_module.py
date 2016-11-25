# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 11:29:56 2016

@author: stefano
"""

#++++++++++++++ Requested Package(s) Import +++++++++++++++#
import matrix_configure_module
import pandas as pd
import numpy as np


#++++++++++++++++++++++ Global Vars +++++++++++++++++++++++#
#verbose = matrix_configure_module.verbose


#++++++++++++++++++++++ Global Funcs ++++++++++++++++++++++#
verbosePrint = matrix_configure_module.verbosePrint
humanSorted = matrix_configure_module.humanSorted


#+++++++++++++++++++++++++++++++++++++++ FUNCTIONS ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

################################################################################
### BUILD MATRIX-LIKE DATAFRAMES (ISs as rows, SAMPLES as columns) FROM ANY_DF #
################################################################################


def buildSeqCountMatrix(any_df, sample_column='association_ID'):
    # Check for required columns
    required_columns = set([sample_column, 'genomic_coordinates', 'seq_count'])
    if required_columns.intersection(set(any_df.columns)) != required_columns:
        verbosePrint("[Warning] Columns {required_columns} are required for SeqCountMatrix. Given data have {columns_found}.".format(required_columns=str(humanSorted(list(required_columns))), columns_found=str(humanSorted(list(any_df)))))
        verbosePrint("...SKIP!")
        return None
    # pivot to get SEQUENCE COUNTS of ISs in df
    seqCount_matrix = pd.pivot_table(any_df, index='genomic_coordinates', columns=sample_column, values='seq_count', aggfunc='sum', margins=False)
    return seqCount_matrix


def buildShsCountMatrix(any_df, sample_column='association_ID'):
    # Check for required columns
    required_columns = set([sample_column, 'genomic_coordinates', 'shearsite'])
    if required_columns.intersection(set(any_df.columns)) != required_columns:
        verbosePrint("[Warning] Columns {required_columns} are required for ShsCountMatrix. Given data have {columns_found}.".format(required_columns=str(humanSorted(list(required_columns))), columns_found=str(humanSorted(list(any_df)))))
        verbosePrint("...SKIP!")
        return None
    # pivot to get counts of distinct SHEARSITES for each ISs in df (unique method exploited by lambda in aggfunc is a method of pd.Series class)
    # NOTE: here MARGINS are not as expected (not partial 'sums' of course, but the 'overall' len(x.unique()) column!!)
    ShsCount_matrix = pd.pivot_table(any_df, index='genomic_coordinates', columns=sample_column, values='shearsite', aggfunc=lambda x: len(x.unique()), margins=False)
    return ShsCount_matrix

    
def buildBarcodeCountMatrix(any_df, sample_column='association_ID'):
    # try to take required columns -> new DF
    required_columns = [sample_column, 'genomic_coordinates', 'randomBC']
    l = []
    try:
        [l.append(any_df.loc[:,c].to_frame()) for c in required_columns]
    except:
        verbosePrint("[Warning] Columns {required_columns} are required for BarcodeCountMatrix. Given data have {columns_found}.".format(required_columns=str(humanSorted(required_columns)), columns_found=str(humanSorted(list(any_df)))))
        verbosePrint("...SKIP!")
        return None
    DF = pd.concat(l, axis=1, join='inner')
    # group-by remaining fields but randomBC - > duplicate entries may arise
    grouped = DF.groupby([sample_column,'genomic_coordinates'], as_index=False)
    # aggregate duplicates taking the 'count-distinct' as values
    barcodeCount_df = grouped.aggregate(lambda x: len(np.unique(x)))
    barcodeCount_df.rename(columns={'randomBC': 'randomBC_countDistinct'}, inplace=True)
    # pivot to shape barcodeCount_df as matrix of distinct RANDOM-BARCODES for each ISs in df
    barcodeCount_matrix = pd.pivot_table(barcodeCount_df, index='genomic_coordinates', columns=sample_column, values='randomBC_countDistinct')
    return barcodeCount_matrix

  
def buildCellCountMatrix(any_df, sample_column='association_ID'):
    # Check for required columns
    required_columns = set([sample_column, 'genomic_coordinates', 'shearsite', 'randomBC'])  # 'shearsite' is redundant, just to be sure about the completeness of the hierarchy (exploited in pivot_table)
    if required_columns.intersection(set(any_df.columns)) != required_columns:
        verbosePrint("[Warning] Columns {required_columns} are required for CellCountMatrix. Given data have {columns_found}.".format(required_columns=str(humanSorted(list(required_columns))), columns_found=str(humanSorted(list(any_df)))))
        verbosePrint("...SKIP!")
        return None
    # pivot to get counts of distinct SHEARSITE-RANDOMTAG couples for each ISs in df
    cellCount_matrix = pd.pivot_table(any_df, index='genomic_coordinates', columns=sample_column, values='randomBC', aggfunc=lambda x: len(x), margins=False)
    return cellCount_matrix


def buildFragmentEstimateMatrix(any_df, sample_column='association_ID'):
    # import rpy2 and sonicLength package
    try:
        import rpy2.robjects as robjects
    except Exception, err_message:
        verbosePrint("[Warning] buildFragmentEstimateMatrix can't run due to some troubles with 'import rpy2.robjects as robjects'.")
        verbosePrint("'import rpy2.robjects as robjects' returned: "+str(err_message))
        verbosePrint("...SKIP!")
        return None
    try:
        from rpy2.robjects.packages import importr
    except Exception, err_message:
        verbosePrint("[Warning] buildFragmentEstimateMatrix can't run due to some troubles with 'from rpy2.robjects.packages import importr'.")
        verbosePrint("'from rpy2.robjects.packages import importr' returned: "+str(err_message))
        verbosePrint("...SKIP!")
        return None
    try:
        sonicLength = importr("sonicLength")
    except Exception, err_message:
        verbosePrint("[Warning] buildFragmentEstimateMatrix can't run due to some troubles occurred while loading 'sonicLength' R-library.")
        verbosePrint('''importr("sonicLength") returned: '''+str(err_message))
        verbosePrint("...SKIP!")
        return None
    # try to take required columns -> new DF
    required_columns = [sample_column, 'genomic_coordinates', 'shearsite']
    l = []
    try:
        [l.append(any_df.loc[:,c].to_frame()) for c in required_columns]
    except:
        verbosePrint("[Warning] Columns {required_columns} are required for FragmentEstimateMatrix. Given data have {columns_found}.".format(required_columns=str(humanSorted(required_columns)), columns_found=str(humanSorted(list(any_df)))))
        verbosePrint("...SKIP!")
        return None
    DF = pd.concat(l, axis=1, join='inner')
    # detect samples
    sample_list = humanSorted(list(DF[sample_column].unique()))
    # list of dataframes with results
    results_list = []
    # loop over samples
    for sample in sample_list:
        try:
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
        except Exception, err_message:
            verbosePrint("  [Warning] Troubles have occurred while running sonicLength R-package.")
            verbosePrint("  error details: "+str(err_message))
            verbosePrint("  ...ABORTED!")
            return None
    # concat dataframes in results_list
    results_df = pd.concat(results_list)
    # pivot to shape results_df as matrix of ESTIMATES OF N PARENTAL FRAGMENT YIELDED BY SONICLENGTH for each ISs in df
    fragmentEstimate_matrix = pd.pivot_table(results_df, index='genomic_coordinates', columns=sample_column, values='fragmentEstimate_count')
    return fragmentEstimate_matrix

