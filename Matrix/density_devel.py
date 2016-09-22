# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 09:41:52 2016

@author: stefano
"""

#++++++++++++++ Requested Package(s) Import +++++++++++++++#
import sys
import pandas as pd
import matrix_configure_module

#++++++++++++++++++++++ Global Vars +++++++++++++++++++++++#
#verbose = matrix_configure_module.verbose


#++++++++++++++++++++++ Global Funcs ++++++++++++++++++++++#
verbosePrint = matrix_configure_module.verbosePrint
humanSorted = matrix_configure_module.humanSorted




def get_ensembles (any_df, per_sample=False, max_dist=3, max_span=5):
    # if per_sample is True, esnambles are computed sample-by-sample, not overall (False)
    # if (locus[n+1] - locus[n] > max_dist), the ensamble is truncated at locus[n]
    # if (locus[n+1] - locus[0] > max_span), the ensamble is truncated at locus[n]
    # NOTE/TrickyUsages: > max_dist=0 OR max_span=1 yield single-covered-base ensambles
    #                    > max_span=0 yield single-ENTRY ensambles (usually NON-SENSE)
    
    verbosePrint("    > Building covered base ensembles ...")
    
    # Add and cast 'chr', 'locus' and 'strand' columns
    verbosePrint("      * Preparing data ...")
    tmp_df = pd.concat([any_df, any_df['genomic_coordinates'].str.split('_').apply(pd.Series)], axis=1)
    tmp_df.rename(columns={0:'chr', 1:'locus', 2:'strand'}, inplace=True)
    tmp_df['locus'] = tmp_df['locus'].astype(int)
    verbosePrint("        Done. {n} total entries.".format(n=str(len(tmp_df))))
    
    # Build groups
    grouping_rule = ['chr', 'strand']
    if per_sample is True:
        grouping_rule += ['association_ID']
    verbosePrint("      * Grouping data by {grouping_rule} ...".format(grouping_rule=str(grouping_rule)))
    grouped_tmp_df = tmp_df.groupby(grouping_rule, as_index=False)
    verbosePrint("        Done. {m} groups.".format(m=str(len(grouped_tmp_df))))
    
    # Iterate over groups: build ensembles
    ensembles = []
    # key[0] -> grouping_rule[0] 'chr'
    # key[1] -> grouping_rule[1] 'strand'
    verbosePrint("      * Looping over groups ...")
    i = 1
    for key, frame in grouped_tmp_df:
        
        # User print and loop count
        verbosePrint("        > {i} of {m} - group {key} - {n} entries to process ...".format(key=str(key), i=str(i), m=str(len(grouped_tmp_df)), n=str(len(frame))))
        i += 1
        
        # Sort/compute-distance strand-wise: add 'dist' column to frame
        if key[1] == '+':
            #verbosePrint("          sorting ...")
            frame = frame.sort(columns='locus', axis=0, ascending=True, inplace=False, kind='quicksort')
            #verbosePrint("          computing ...")
            frame['dist'] = frame.loc[:,'locus'].shift(-1) - frame.loc[:,'locus']
        elif key[1] == '-':
            #verbosePrint("          sorting ...")
            frame = frame.sort(columns='locus', axis=0, ascending=False, inplace=False, kind='quicksort')
            #verbosePrint("          computing ...")
            frame['dist'] = frame.loc[:,'locus'] - frame.loc[:,'locus'].shift(-1)
        else:
            print "\n[ERROR] Unsupported 'strand' symbol found: '{strand}'.".format(strand=str(key[1]))
            sys.exit("\n[QUIT]\n")
        
        # Iterate over frame's entries to build ensambles
        frame_as_ensambles = []
        frame_as_ensambles_append = frame_as_ensambles.append
        frame.reset_index(inplace=True, drop=True)  # mandatory to exploit 'index' in loop with 'iloc'
        last_index = 0
        span = 1
        #verbosePrint("          splitting ...")
        for index, row in frame.iterrows():
            
            # mandatory to close the loop appending last ensamble, and needed to handle 1-entry limit-case
            if pd.isnull(row['dist']):
                ens = frame.iloc[last_index:index+1,:]
                ens = ens.drop(['chr', 'locus', 'strand', 'dist'], axis=1)
                frame_as_ensambles_append(ens)
            # extract chunk (ensamble) and append to frame_as_ensambles list
            elif (row['dist'] > max_dist) or (span + row['dist'] > max_span):
                ens = frame.iloc[last_index:index+1,:]
                ens = ens.drop(['chr', 'locus', 'strand', 'dist'], axis=1)
                frame_as_ensambles_append(ens)
                last_index = index + 1
                span = 1
            # just update span and next!
            else:
                span += row['dist']
        
        # join local frame_as_ensambles list to global ensambles list
        #verbosePrint("          appending results ...")
        ensembles += frame_as_ensambles
        verbosePrint("          {z} ensembles built.".format(z=str(len(frame_as_ensambles))))
    
    verbosePrint("    > Done! {tot} ensembles built.".format(tot=str(len(ensembles))))
    return ensembles




def classic_IS (ensemble, place_on_mode=True, **kwargs):
    # NOTE: Kwargs are needed only to verbosePrint operations: if provided, must be both "n" and "tot" (in the spirit of "processing n out of tot ...")
    
    if kwargs:
            verbosePrint("        > processing {n} of {tot} ... ".format(n=str(kwargs['n']), tot=str(kwargs['tot'])))
            
    # Fast exit for single-base ensembles
    if len(set(ensemble['genomic_coordinates'])) < 2:
        verbosePrint("          Fast IS: {is_id}".format(is_id=str(list(set(ensemble['genomic_coordinates']))[0])))
        return ensemble.copy()
    verbosePrint("          Slow IS: {is_ids}".format(is_ids=str(set(ensemble['genomic_coordinates']))))
    
    # Get coordinates and seq_count in tmp_df
    tmp_df = pd.concat([ensemble['seq_count'], ensemble['genomic_coordinates'].str.split('_').apply(pd.Series).rename(columns={0:'chr', 1:'locus', 2:'strand'})], axis=1)
    tmp_df['locus'] = tmp_df['locus'].astype(int)
    
    # Find IS_locus
    IS_locus = None
    if place_on_mode is True:
        IS_locus = tmp_df.loc[tmp_df['seq_count'].idxmax(), 'locus']  # seq_count max
    else:
        IS_locus = tmp_df.iloc[0]['locus']  # the first (tmp_df is already properly sorted - via ensemble)
    
    # Create IS_genomic_coordinate
    IS_genomic_coordinate = "_".join([tmp_df.iloc[0]['chr'], str(IS_locus), tmp_df.iloc[0]['strand']])
    
    # Compute shearsite_correction (strand-wise, always sum to correct!), if shearsite is present. 
    if 'shearsite' in ensemble.columns:
        if tmp_df.iloc[0]['strand'] == '+':
            tmp_df['shearsite_diff'] = tmp_df['locus'] - IS_locus
        elif tmp_df.iloc[0]['strand'] == '-':
            tmp_df['shearsite_diff'] = -1 * (tmp_df['locus'] - IS_locus)
        else:
            print "\n[ERROR] Unsupported 'strand' symbol found: '{strand}'.".format(strand=str(tmp_df.iloc[0]['strand']))
            sys.exit("\n[QUIT]\n")
    
    # instantiate raw_IS
    raw_IS = ensemble.copy()
    
    # fix 'genomic_coordinates'
    raw_IS['genomic_coordinates'] = IS_genomic_coordinate
    
    # correct 'shearsite' if needed
    if 'shearsite_diff' in tmp_df.columns:
        raw_IS['shearsite'] = raw_IS['shearsite'] + tmp_df['shearsite_diff']
    
    # group-by and aggregate to fix the mess -> IS
    grouping_rule = [x for x in raw_IS.columns if x not in ['header_list', 'seq_count']]
    grouped = raw_IS.groupby(grouping_rule, as_index=False)
    
    # aggregate
    IS = grouped.agg('sum')  ### BUGGY BECAUSE 'header_list' col is not returned aggregated by sum (list concatenation)
    
    return IS




def compute_ISs (any_df, method='classic'):
    
    # Compute ensembles
    ensembles = get_ensembles (any_df)
    
    # Get ISs (list of any_df chunks properly processed)
    verbosePrint("    > Computing ISs over covered base ensembles ...")
    ISs = None
    if method=='classic':
        verbosePrint("      * method: 'classic'")
        if 'header_list' in any_df.columns:
            verbosePrint("        [WARNING] 'classic' method for ISs computing is still buggy and may also CRUSH while trying handle 'header_list'!")
        #ISs = [classic_IS(e) for e in ensembles]
        ISs = [classic_IS(e, n=n, tot=len(ensembles)) for n,e in enumerate(ensembles, start=1)]
    else:
        print "\n[ERROR] Unsupported 'method={m}' for compute_ISs".format(m=str(method))
        sys.exit("\n[QUIT]\n")
    verbosePrint("    > Done! {tot} ISs built.".format(tot=str(len(ISs))))
    
    # Reconstruct data
    verbosePrint("    > Reshaping data ...")
    ISs_any_df = pd.concat(ISs)
    ISs_any_df.reset_index(inplace=True, drop=True)
    verbosePrint("      Done! DataFrame built. {n} entries.".format(n=str(len(ISs_any_df))))
    
    return ISs_any_df
            
    











    