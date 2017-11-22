# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 10:15:22 2017

@author: Adriano De Marino
"""

#++++++++++++++++++++++++++++ Requested Package(s) Import +++++++++++++++++++++++++++++#

from datetime import date
import itertools
import pandas as pd
import numpy as np
from collections import Counter

#+++++++++++++++++++++++++++++++++++++++ FUNCTIONS +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

def pivot(any_df):
    '''
    The function pivot return a matrix like this:
                                             sample1    sample2    sample3  ...
    chr locus strand shearsite randomBC 
    chr1  123551  +   223     GAAATTGCCCGT      1                   1
    ....
    the information from chr to randomBC are the essential of the concept of UMI (unique molecular identifier)
    ....
    '''
    new=any_df
    new[['chr','locus','strand']] = (new.genomic_coordinates.str.split('_', expand= True ))
    matrix_UMI = new.pivot_table(index=['chr','locus','strand','shearsite','randomBC'], columns=['association_ID'],values='seq_count',fill_value=0)
    return matrix_UMI

def contamina(dataframe_UMI,metadata):
    
    _stendino = []
    for index, columns in dataframe_UMI.iterrows():
        Who_have_it = np.flatnonzero(columns)
        if len(Who_have_it) > 1:
            #combinazioni a coppie in un array exp. ([1,2,55]) --> [[1,2],[1,55],[2,55]]
            combinations = list(itertools.combinations(Who_have_it, 2))
            #ciclo sulle sull'array di combinazioni per realizzare confronti a coppie.
            for pair in combinations:
                _sample_name_1, _sample_name_2 = dataframe_UMI.columns[pair[0]], dataframe_UMI.columns[pair[1]]
                _seq_count_1, _seq_count_2 = columns[pair[0]], columns[pair[1]]
                _info_1, _info_2 = metadata[metadata['CompleteAmplificationID'] == _sample_name_1], metadata[metadata['CompleteAmplificationID'] == _sample_name_2]
                _dateOne = map(int, _info_1[['FusionPrimerPCRDate']].values[0][0].split('-'))
                _date2 = map(int, _info_2[['FusionPrimerPCRDate']].values[0][0].split('-'))
                _LTR_LINKER_1, _LTR_LINKER_2 = _info_1[['TagID']].values[0][0], _info_2[['TagID']].values[0][0]
                _LTR_1, _LTR_2 = _LTR_LINKER_1.split('.')[0], _LTR_LINKER_2.split('.')[0]
                _LINKER_1, _LINKER_2 = _LTR_LINKER_1.split('.')[1], _LTR_LINKER_2.split('.')[1]
                x0 = date(*_dateOne)
                x1 = date(*_date2)
                _distance = x0 - x1
                
                if _LTR_1 != _LTR_2 and _LINKER_1 == _LINKER_2:
                    Type = 'VP'
                    Contamination = 'YES'
                elif [ _LTR_1 != _LTR_2 or _LTR_1 == _LTR_2 ] and _LINKER_1 != _LINKER_2:
                    Type = 'FP'
                    Contamination = 'NO'
                else:
                    Type = 'NONE'
                    Contamination = 'YES'

                adjast = []
                adjast.extend(index)
                adjast.extend(( _sample_name_1, _seq_count_1, _sample_name_2, _seq_count_2, _distance.days, _LTR_LINKER_1, _LTR_LINKER_2, Type, Contamination))
                _stendino.append(adjast)
    
    return pd.DataFrame(_stendino, columns=('chr', 'locus', 'strand', 'shearsite', 'randomBC', 'X', 'seqCount_X', 'Y', 'seqCount_Y', 'DayDistance', 'TagID_X', 'TagID_Y', 'Type', 'True_Contamination?' ))

def comparison(dataframe_UMI):
    '''This function return a numpy.matrix in which columns and rows there are pooled samples
            in each cells of matrix there are the number of umi shared from one sample with the other 
            in the diagonal instead there is the number of UMI own of a specific sample.'''
    matrix = []
    for column in dataframe_UMI:
        #retrieve the UMI > 0 in current sample
        current_sample = dataframe_UMI.loc[ dataframe_UMI[column] > 0 ]
        #sum the number of no empty value
        current_sample = (current_sample != 0).astype(int).sum()
        #creation of array where the 1 reppresent the current sample, the value not blank instead the number of ratio between shared element / element of the current sample 
        matrix.append(current_sample.values)
    del current_sample
    return matrix

def skip_diag_strided(A):
    '''This function return a matrix without diagonal gived in input a matrix'''
    m = A.shape[0]
    strided = np.lib.stride_tricks.as_strided
    s0,s1 = A.strides
    return strided(A.ravel()[1:], shape=(m-1,m), strides=(s0+s1,s1)).reshape(m,-1)

def isSymmetric(matrix, N):
    for i in range(N):
        for j in range(N):
            if (matrix[i][j] != matrix[j][i]):
                return False
    return True

def associationFILE(file_metadata, ProjectID):
    metadata = pd.read_csv(file_metadata, sep='\t')
    # listP = ProjectID.strip().replace(' ','').split(',')
    for i in ProjectID:
        arg = metadata[metadata['ProjectID'] == i]
        metadata = metadata.append(arg)
    return metadata[['TagID','FusionPrimerPCRDate','CompleteAmplificationID']]

def describe(matrix,df):
    '''This function return a SUMMARY information about a square simmetric Matrix'''
    N=len(matrix)
    'Symmetric Matrix' if isSymmetric(matrix,N) == True else exit()

    #triangolar upper matrix with diagonal
    triu = np.triu(matrix,0)
    #ratio of triangolar upper matrix with diagonal
    triu_ratio_YES_diagonal = np.true_divide(triu,triu.diagonal()[:,None])
    #ratio of triangolar upper matrix NO diagonal
    triu_ratio_NO_diagonal = skip_diag_strided(triu_ratio_YES_diagonal)

    #list of shared mean --> each element into array is the average between one sample with the others
    sharingMedia_byRow = triu_ratio_NO_diagonal.mean(1)

    sharingMedia_noZero = sharingMedia_byRow[np.nonzero(sharingMedia_byRow)] #solo medie dei campioni realmente contaminati

    contaminated_sample = len(sharingMedia_noZero)
    frequency = Counter(sharingMedia_noZero)
    list_contaminated_sample = [] ; [ list_contaminated_sample.append(df.columns[i]) if sharingMedia_byRow[i] != 0 else '' for i in range(0,len(sharingMedia_byRow))]
    
    if len(sharingMedia_noZero) != 0:
        Contamination_rate_LOWER = np.min(sharingMedia_noZero)
        Frequency_of_contamination_rate_LOWER = frequency[np.min(sharingMedia_noZero)]
        The_contaminated_SAMPLE_shared_the_MAIN_number_of_UMI = df.columns[sharingMedia_byRow.argmax()]
        The_contaminated_SAMPLE_shared_the_smallest_number_of_UMI = df.columns[np.where(sharingMedia_byRow == sharingMedia_noZero.min())[0][0]]
        Average_share_of_only_contaminated_samples = sharingMedia_noZero.mean()
        LIST_of_contaminated_SAMPLES = list_contaminated_sample
        median_only_contaminated = np.median(sharingMedia_noZero)
        percentile = np.percentile(sharingMedia_noZero,[25,50,75])
    else:
        Contamination_rate_LOWER = 0
        Frequency_of_contamination_rate_LOWER = 0
        The_contaminated_SAMPLE_shared_the_MAIN_number_of_UMI = 'NONE'
        The_contaminated_SAMPLE_shared_the_smallest_number_of_UMI = 'NONE'
        LIST_of_contaminated_SAMPLES = 'NONE'
        Average_share_of_only_contaminated_samples = 0
        median_only_contaminated = 'NONE'
        percentile = 'NONE'


    info_list = [ "Number_of_TOTAL_UMI","Number_of_SHARED_UMI","Ratio_of_SHARED/TOTAL_UMI","Number_of_TOTAL_contaminated_Samples","Contamination_rate_HIGHER", \
    			  "Frequency_of_contamination_rate_HIGHER","Contamination_rate_LOWER","Frequency_of_contamination_rate_LOWER","LIST_of_contaminated_SAMPLES", \
    			  "The_contaminated_SAMPLE_shared_the_MAIN_number_of_UMI","The_contaminated_SAMPLE_shared_the_smallest_number_of_UMI","Average_share_of_all_samples", \
    			  "Average_share_of_only_contaminated_samples","median_all_sample","median_only_contaminated","percentile"
    			]

    Number_of_TOTAL_UMI                                       = np.trace(triu)
    Number_of_SHARED_UMI                                      = np.triu(matrix,1).sum()
    Ratio_of_SHARED_TOTAL_UMI                                 = Number_of_SHARED_UMI/float(Number_of_TOTAL_UMI) * 100
    Number_of_TOTAL_contaminated_Samples                      = contaminated_sample
    Contamination_rate_HIGHER                                 = np.max(sharingMedia_byRow)
    Frequency_of_contamination_rate_HIGHER                    = frequency[np.max(sharingMedia_byRow)]
    Average_share_of_all_samples                              = triu_ratio_NO_diagonal.mean()
    median_all_sample                                         = np.median(sharingMedia_byRow)
    quantili                                                  = percentile


    cmd = [ Number_of_TOTAL_UMI, Number_of_SHARED_UMI, Ratio_of_SHARED_TOTAL_UMI, Number_of_TOTAL_contaminated_Samples, Contamination_rate_HIGHER, Frequency_of_contamination_rate_HIGHER, Contamination_rate_LOWER, Frequency_of_contamination_rate_LOWER, LIST_of_contaminated_SAMPLES, The_contaminated_SAMPLE_shared_the_MAIN_number_of_UMI, The_contaminated_SAMPLE_shared_the_smallest_number_of_UMI, Average_share_of_all_samples,Average_share_of_only_contaminated_samples, median_all_sample, median_only_contaminated, quantili]

    s = []

    for i,j in zip(info_list,cmd):
        s.append((i,j))
    return pd.DataFrame(s, columns=('info','data'))
   		






















