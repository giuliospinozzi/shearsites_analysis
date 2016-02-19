# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 16:39:41 2016

@author: stefano
"""

### PACKAGES IMPORT ##############################################
import matrix_RandomBC_outputModule
import matrix_RandomBC_globModule
import os
##################################################################

### IMPORT GLOBAL VARS ###########################################################
# Print
verbose = matrix_RandomBC_globModule.verbose
# Output
common_output_ground_dir = matrix_RandomBC_globModule.common_output_ground_dir
cem_data_outfolder = matrix_RandomBC_globModule.cem_data_outfolder
cem_data_outfile_name = matrix_RandomBC_globModule.cem_data_outfile_name
##################################################################################

### IMPORT GLOBAL FUNCS #########################################
verbosePrint = matrix_RandomBC_globModule.verbosePrint
humanSorted = matrix_RandomBC_globModule.humanSorted
#################################################################

### FIXED LOCAL VARIABLES - CEM CONFIG ###############################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################         
cem_locations = [('11',64537165,64537171), ('17',2032349,2032355), ('17',47732336,47732342), ('2',24546568,24546574), ('2',73762395,73762401), ('16',28497495,28497501), ('3', 52306691, 52306697), ('8', 8866483, 8866489)]  # [('chormosome', min_accepted_locus, max_accepted_locus), ... ]
cem_strands = ['-', '-', '-', '-', '-', '-', '-', '+']
cem_symbols_list = ['CEM6_37_1','CEM6_37_2','CEM6_37_3','CEM6_37_4','CEM6_37_5','CEM6_37_6', 'CEM1_25', 'CEM1_41_illumina']
cem_coordinates_list = ['chr11:64537168:-','chr17:2032352:-','chr17:47732339:-','chr2:24546571:-','chr2:73762398:-','chr16:28497498:-', 'chr3:52306694:-', 'chr8:8866486:+']
cem_quantif_expect_dict = {('CEM1_41_illumina', 'L'): 0.28, ('CEM6_37_1', 'L'): 0.12, ('CEM6_37_2', 'L'): 0.12, ('CEM6_37_3', 'L'): 0.12, ('CEM6_37_4', 'L'): 0.12, ('CEM6_37_5', 'L'): 0.12, ('CEM6_37_6', 'L'): 0.12, ('CEM1_41_illumina', 'M'): 0.24, ('CEM6_37_1', 'M'): 0.12, ('CEM6_37_2', 'M'): 0.12, ('CEM6_37_3', 'M'): 0.12, ('CEM6_37_4', 'M'): 0.12, ('CEM6_37_5', 'M'): 0.12, ('CEM6_37_6', 'M'): 0.12, ('CEM1_41_illumina', 'N'): 0.2, ('CEM6_37_1', 'N'): 0.12, ('CEM6_37_2', 'N'): 0.12, ('CEM6_37_3', 'N'): 0.12, ('CEM6_37_4', 'N'): 0.12, ('CEM6_37_5', 'N'): 0.12, ('CEM6_37_6', 'N'): 0.12, ('CEM1_41_illumina', 'O'): 0.1, ('CEM6_37_1', 'O'): 0.12, ('CEM6_37_2', 'O'): 0.12, ('CEM6_37_3', 'O'): 0.12, ('CEM6_37_4', 'O'): 0.12, ('CEM6_37_5', 'O'): 0.12, ('CEM6_37_6', 'O'): 0.12, ('CEM1_41_illumina', 'P'): 0.05, ('CEM6_37_1', 'P'): 0.12, ('CEM6_37_2', 'P'): 0.12, ('CEM6_37_3', 'P'): 0.12, ('CEM6_37_4', 'P'): 0.12, ('CEM6_37_5', 'P'): 0.12, ('CEM6_37_6', 'P'): 0.12, ('CEM1_41_illumina', 'Q'): 0.03, ('CEM6_37_1', 'Q'): 0.12, ('CEM6_37_2', 'Q'): 0.12, ('CEM6_37_3', 'Q'): 0.12, ('CEM6_37_4', 'Q'): 0.12, ('CEM6_37_5', 'Q'): 0.12, ('CEM6_37_6', 'Q'): 0.12, ('CEM1_41_illumina', 'R'): 0.0029994, ('CEM6_37_1', 'R'): 0.119976, ('CEM6_37_2', 'R'): 0.119976, ('CEM6_37_3', 'R'): 0.119976, ('CEM6_37_4', 'R'): 0.119976, ('CEM6_37_5', 'R'): 0.119976, ('CEM6_37_6', 'R'): 0.119976, ('CEM1_41_illumina', 'S'): 0.00030003, ('CEM6_37_1', 'S'): 0.120012, ('CEM6_37_2', 'S'): 0.120012, ('CEM6_37_3', 'S'): 0.120012, ('CEM6_37_4', 'S'): 0.120012, ('CEM6_37_5', 'S'): 0.120012, ('CEM6_37_6', 'S'): 0.120012, ('CEM1_41_illumina', 'T'): 0, ('CEM6_37_1', 'T'): 0.12, ('CEM6_37_2', 'T'): 0.12, ('CEM6_37_3', 'T'): 0.12, ('CEM6_37_4', 'T'): 0.12, ('CEM6_37_5', 'T'): 0.12, ('CEM6_37_6', 'T'): 0.12, ('CEM1_41_illumina', 'UT'): 0, ('CEM6_37_1', 'UT'): 0, ('CEM6_37_2', 'UT'): 0, ('CEM6_37_3', 'UT'): 0, ('CEM6_37_4', 'UT'): 0, ('CEM6_37_5', 'UT'): 0, ('CEM6_37_6', 'UT'): 0}
######################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

### OTHER LOCAL VARIABLES ####################################################################
# Output
sep = '\t'
eol = '\n'
# misc
relabelling_sep = "@@"  # choose a robust one, it's just exploited in computation
##############################################################################################


#+++++++++++++++++++++++++++++++++++++++ FUNCTIONS ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

def buildIndexes(cem_locations_tuple, strand):
    chrom = cem_locations_tuple[0]
    start_region = cem_locations_tuple[1]
    end_region = cem_locations_tuple[2]
    indexes = []
    for locus in range(start_region, end_region+1):
        indexes.append('_'.join([chrom, str(locus), strand]))
    return indexes
    
def extractCEM_fromMatrix(matrix_df, cem_locations=cem_locations, cem_strands=cem_strands):
    indexes = []
    for i in range(len(cem_locations)):
        indexes += buildIndexes(cem_locations[i], cem_strands[i])
    cem_matrix_df = matrix_df.loc[indexes]
    cem_matrix_df.dropna(axis=0, how='all', inplace=True)
    return cem_matrix_df
    
def extractCEM_fromDF(any_df, cem_locations=cem_locations, cem_strands=cem_strands):
    coordinates = []
    for i in range(len(cem_locations)):
        coordinates += buildIndexes(cem_locations[i], cem_strands[i])
    cem_df = any_df[any_df['genomic_coordinates'].isin(coordinates)]
    return cem_df

def CEM_info(integration_id, cem_locations=cem_locations, cem_strands=cem_strands, cem_coordinates=cem_coordinates_list, cem_symbols = cem_symbols_list):
    n_cem = None
    for i in range(len(cem_locations)):
        cem_i_indexes = buildIndexes(cem_locations[i], cem_strands[i])
        if integration_id in cem_i_indexes:
            n_cem = i
            break
    info_dict = {'coordinate': cem_coordinates[n_cem],
                 'symbol': cem_symbols[n_cem]}
    return info_dict
    
def CEM_expected_quantification(quantification_id_tuple, cem_quantif_expect_dict=cem_quantif_expect_dict):
    return {'expected quantification': cem_quantif_expect_dict[quantification_id_tuple]}
    
    
def buildRelabellingDict(asso_dict, concat, use_fields=(3,2,4,7,8)):
    relabellingDict = {}
    for k, d in asso_dict.items():
        fields = []
        [fields.append(d[n]) for n in use_fields]
        new_label = concat.join(fields)
        relabellingDict[k] = new_label
    return relabellingDict
    
def relabelling(df, asso_dict, concat=relabelling_sep, inplace=False):
    # NOTE: - df col names must be the keys of asso_dict!
    if inplace is True:
        # copy : boolean, default True
        df.rename(columns=buildRelabellingDict(asso_dict, concat), inplace=True)
    else:
        # new object
        return df.rename(columns=buildRelabellingDict(asso_dict, concat), inplace=False)


#+++++++++++++++++++++++++++++++++++++++ MAIN FUNCTION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#



def exportCEM(df_dict, asso_dict, filename=cem_data_outfile_name, sep=sep, eol=eol):
    
    # e.g. df_dict = {'sequence_count': df_seqCount_matrix,
    #                 'shearsite_count': df_ShsCount_matrix,
    #                 'shearsite&TAG_count': df_cellCount_matrix,
    #                 ...,
    #                }
    # NOTE: matrixes in df_dict can be any, but for a proper comparison with
    # th_abundance they are expected to be normalized!
    
    def arrangeData(df_dict, asso_dict, relabelling_sep=relabelling_sep):
        # e.g. df_dict = {'sequence_count': df_seqCount_matrix,
        #                 'shearsite_count': df_ShsCount_matrix,
        #                 'shearsite&TAG_count': df_cellCount_matrix,
        #                 ...,
        #                }
        # NOTE: matrixes in df_dict can be any, but for a proper comparison with
        # th_abundance they are expected to be normalized!
        export_dict = {}
        for name in humanSorted(df_dict.keys()):
            df = df_dict[name]
            df = extractCEM_fromMatrix(df)
            df = relabelling(df, asso_dict)
            df = df.unstack()
            for index, value in df.iteritems():
                # export_dict structure
                if index not in export_dict.keys():
                    export_dict[index] = {name:value}
                else:
                    export_dict[index][name] = value
                # add cem info (redundant update)
                info_dict = CEM_info(index[-1])
                export_dict[index].update(info_dict)
                # add expected quantification  (redundant update)
                condition_tuple = tuple(index[0].split(relabelling_sep))
                quantification_id_tuple = tuple([info_dict['symbol'], condition_tuple[2]])
                expected_quantification_dict = CEM_expected_quantification(quantification_id_tuple)
                export_dict[index].update(expected_quantification_dict)
        # e.g. export_dict = {('column field concatenation with relabelling_sep', 'integration id / index'): another_dict}    
        #                                       another_dict = {df_dict key1: value, df_dict key2: value, ..., 'coordinate': string, 'symbol': string, 'expected quantification': value}
        return export_dict
        
    def formatData(export_dict, df_dict, relabelling_sep=relabelling_sep, header=True):
        # e.g. export_dict = {('column field concatenation with relabelling_sep', 'integration id / index'): another_dict}    
        #                                       another_dict = {df_dict key1: value, df_dict key2: value, ..., 'coordinate': string, 'symbol': string, 'expected quantification': value}
        lines = []
        quantification_labels = []
        quantification_keys = []
        for quantification_name in humanSorted(df_dict.keys()):
            quantification_labels.append(quantification_name)
            quantification_keys.append(quantification_name)
        if header:
            lines = [['label', 'approach', 'wet_method', 'dilution', 'bio_replicate', 'tech_replicate', 'cem_region', 'cem_genomic_id']+quantification_labels+['th_abundance']]
        for k, d in export_dict.items():
            #index = k[-1].split(relabelling_sep)
            condition_tuple = tuple(k[0].split(relabelling_sep))
            # Get data
            approach = condition_tuple[0]
            wet_method = condition_tuple[1]
            dilution = condition_tuple[2]
            bio_replicate = condition_tuple[3]
            tech_replicate = condition_tuple[4]
            cem_region = d['symbol']
            cem_genomic_id = d['coordinate']
            th_abundance = str(d['expected quantification'])
            label = '_'.join([approach, wet_method, dilution, bio_replicate, tech_replicate, cem_genomic_id])
            quantification_list = []
            for l in quantification_keys:
                quantification_list.append(str(d[l]))
            # create line as list and append to lines
            line = [label, approach, wet_method, dilution, bio_replicate, tech_replicate, cem_region, cem_genomic_id]+quantification_list+[th_abundance]
            lines.append(line)
            # lines = file lines as nested list
        return lines
    
    # Code
    verbosePrint("\nExport CEM data ...")
    export_dict = arrangeData(df_dict, asso_dict)
    lines = formatData(export_dict, df_dict)
    OUTDIR = matrix_RandomBC_outputModule.buildOutputPath(common_output_ground_dir, cem_data_outfolder)
    outpath = os.path.normpath(os.path.join(OUTDIR, filename))
    with open(outpath, 'w') as out_stream:
        rows = [sep.join(l) for l in lines]
        out_stream.write(eol.join(rows)+eol)
    verbosePrint(">>> Created file '{outpath}'".format(outpath=str(outpath)))
    
    return outpath
