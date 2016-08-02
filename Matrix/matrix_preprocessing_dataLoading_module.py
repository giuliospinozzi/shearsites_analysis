# -*- coding: utf-8 -*-
"""
Created on Tue Aug  2 10:15:10 2016

@author: stefano
"""


#++++++++++++++++++++++++++++ Requested Package(s) Import +++++++++++++++++++++++++++++#
import pandas as pd
import matrix_configure_module

from matrix_preprocessing_dataSources_module import getLaunchPathDict
from matrix_preprocessing_dataGathering_module import buildDataFrame
from matrix_preprocessing_dataGathering_module import loadRefactored_asDataframe
from matrix_preprocessing_dataGathering_module import loadFasta_asDataframe


#++++++++++++++++++++++ Global Vars +++++++++++++++++++++++#
verbose = matrix_configure_module.verbose


#++++++++++++++++++++++ Global Funcs ++++++++++++++++++++++#
verbosePrint = matrix_configure_module.verbosePrint
humanSorted = matrix_configure_module.humanSorted


#+++++++++++++++++++++++++++++++++++++++ FUNCTIONS +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#


def loadData (launch_path_dict, drop_headers, compression):
    
    # Check if data are all random barcoded
    def rBC_globally_available (launch_path_dict):
        for DISEASE in launch_path_dict.keys():
            for PATIENT in launch_path_dict[DISEASE].keys():
                for POOL in launch_path_dict[DISEASE][PATIENT].keys():
                    if launch_path_dict[DISEASE][PATIENT][POOL].get('randomBC_data_path') is None:
                        return False
        return True

    def loadPool (pool_dict, drop_headers, force_drop_rBC, compression):
        source_data_path = pool_dict['source_data_path']  # Key error if not present
        refactored_DF = loadRefactored_asDataframe (source_data_path, compression=compression)
        randomBC_data_path = pool_dict.get('randomBC_data_path') # Return None instead of Key error if not present
        # Force to ignore randomBC, even if available
        if force_drop_rBC:
            randomBC_data_path = None
        rBC_fasta_DF = loadFasta_asDataframe (randomBC_data_path, compression=compression)
        any_df = buildDataFrame (refactored_DF, rBC_fasta_DF=rBC_fasta_DF, drop_headers=drop_headers)
        return any_df
    
    verbosePrint("\n>>> Loading data ...")
    # Drop randomBC, unless they are available for all data
    force_drop_rBC = True
    if rBC_globally_available(launch_path_dict):
        force_drop_rBC = False
    else:
        verbosePrint("[Warning] Random Barcodes will be ignored to harmonize data contents")
    # Loop over launch_path_dict and load data pool-by-pool
    any_df_per_pool = []
    for DISEASE in humanSorted(launch_path_dict.keys()):
        verbosePrint("> DISEASE:", DISEASE)
        for PATIENT in humanSorted(launch_path_dict[DISEASE].keys()):
            verbosePrint("\t> PATIENT:", PATIENT)
            for POOL in humanSorted(launch_path_dict[DISEASE][PATIENT].keys()):
                verbosePrint("\t\t> POOL: {POOL}. Loading ...".format(POOL=str(POOL)))
                any_df_per_pool.append(loadPool(launch_path_dict[DISEASE][PATIENT][POOL], drop_headers, force_drop_rBC, compression))
    verbosePrint(">>> Data loaded!")
    # append dataframes
    verbosePrint("\n>>> Building Dataframe ...")
    any_df = pd.concat(any_df_per_pool, axis=0, join='inner', ignore_index=True)
    verbosePrint(">>> Done!")
    # return a unique any_df for all the data
    return any_df


#++++++++++++++++++++++++++++++++++++++ MAIN and TEST +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#


if __name__ == "__main__":
    
    ### Test vars ###
    # 2 dataset:
    #   - Ferrari, thalassemia, all pools
    #   - HIV, circulatingDNA, pool 'cLR1' only
    dataset_tuple_list = [ ('/opt/NGS/results', 'Ferrari', 'PE_Thalassemia_inVitro', False), ('/opt/NGS/results', 'HIV_patients', 'circDNA', ['cLR1']) ]
    drop_headers = True
    compression = 'gzip'
    
    ### Test code ###
    verbosePrint("\n[START] \n")
    verbosePrint("[INPUT VAR]")
    verbosePrint("> dataset_tuple_list: {dataset_tuple_list}".format(dataset_tuple_list=str(dataset_tuple_list)))
    verbosePrint("> drop_headers: {drop_headers}".format(drop_headers=str(drop_headers)))
    verbosePrint("> compression: {compression}".format(compression=str(compression)))
    # Load data Paths
    verbosePrint("\n[LOADING PATHS]")
    launch_path_dict = getLaunchPathDict(dataset_tuple_list)
    # Load Data
    verbosePrint("\n[LOADING DATA]")
    any_df = loadData (launch_path_dict, drop_headers, compression)
    verbosePrint("\n[OUTPUT FEATURES]")
    verbosePrint("> len: {len_df}".format(len_df=str(len(any_df))))
    verbosePrint("> columns: {col}".format(col=str(list(any_df.columns))))
    verbosePrint("\n[OUTPUT HEAD]")
    verbosePrint(any_df.head(20))
    verbosePrint("\n[OUTPUT TAIL]")
    verbosePrint(any_df.tail(20))
    
#    # Write data to check them
#    import os
#    out_path = os.path.join(os.getcwd(), "data_dump.csv")
#    verbosePrint("\n>>> Export Dataframe to check data (path = '{out_path}' ) ...".format(out_path=str(out_path)))
#    any_df.to_csv(path_or_buf=out_path, sep='\t', encoding='utf-8')
#    verbosePrint( ">>> File Created!"
    
    verbosePrint("\n[QUIT] \n")
    
