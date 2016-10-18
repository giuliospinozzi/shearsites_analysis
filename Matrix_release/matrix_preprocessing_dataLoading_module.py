# -*- coding: utf-8 -*-
"""
Created on Tue Aug  2 10:15:10 2016

@author: stefano
"""


#++++++++++++++++++++++++++++ Requested Package(s) Import +++++++++++++++++++++++++++++#
import pandas as pd
import matrix_configure_module

from matrix_preprocessing_dataGathering_module import buildDataFrame
from matrix_preprocessing_dataGathering_module import loadRefactored_asDataframe
from matrix_preprocessing_dataGathering_module import loadFasta_asDataframe


#++++++++++++++++++++++ Global Vars +++++++++++++++++++++++#
#verbose = matrix_configure_module.verbose


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
        verbosePrint("      > Loading refactored ...")
        refactored_DF = loadRefactored_asDataframe (source_data_path, compression=compression)
        verbosePrint("        done!")
        randomBC_data_path = pool_dict.get('randomBC_data_path') # Return None instead of Key error if not present
        # Force to ignore randomBC, even if available
        if force_drop_rBC:
            randomBC_data_path = None
        if randomBC_data_path is not None:
            verbosePrint("      > Loading random barcodes ...")
        rBC_fasta_DF = loadFasta_asDataframe (randomBC_data_path, compression=compression)
        if randomBC_data_path is not None:
            verbosePrint("        done!")
        verbosePrint("      > Shaping data as dataframe ...")
        any_df = buildDataFrame (refactored_DF, rBC_fasta_DF=rBC_fasta_DF, drop_headers=drop_headers)
        verbosePrint("        done!")
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
        verbosePrint("> DISEASE: '{DISEASE}'".format(DISEASE=str(DISEASE)))
        for PATIENT in humanSorted(launch_path_dict[DISEASE].keys()):
            verbosePrint("  > PATIENT: '{PATIENT}'".format(PATIENT=str(PATIENT)))
            for POOL in humanSorted(launch_path_dict[DISEASE][PATIENT].keys()):
                verbosePrint("    > POOL: '{POOL}'".format(POOL=str(POOL)))
                any_df_per_pool.append(loadPool(launch_path_dict[DISEASE][PATIENT][POOL], drop_headers, force_drop_rBC, compression))
    verbosePrint(">>> Data loaded!")
    # append dataframes
    verbosePrint("\n>>> Joining Dataframe(s) ...")
    any_df = pd.concat(any_df_per_pool, axis=0, join='inner', ignore_index=True)
    verbosePrint(">>> Done!")
    # return a unique any_df for all the data
    return any_df

