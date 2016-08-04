# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 15:16:12 2016

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


#+++++++++++++++++++++++++++++++++++++++ FUNCTIONS ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

### NOTE ABOUT findEDthreshold(N):
### May be exploited by filterBy_randomBC_EditDistance as kwarg, e.g. "ED_rule = findEDthreshold";
### in this case, the code of filterBy_randomBC_EditDistance should be changed from
### 'ACTUAL VERSION - CODE 1' TO 'EXPERIMENTAL VERSION - CODE 2'

#def findEDthreshold(N):
#    ### INPUT:
#    # N is the number of distinct randomBarcodes (usually within a given ShS)
#    # p[x] is the probability that 2 randomBarcodes have ED = x
#    EDthreshold = 5  # 9, the expected value, -3 (variance rounded in excess), -1 (just by design of filtering algorithm)
#    p = [5.96e-09, 2.145767e-06, 3.540516e-05, 0.0003540516, 0.002389848, 0.01147127, 0.04014945, 0.1032414, 0.1935777, 0.2581036, 0.2322932, 0.1267054, 0.03167635]
#    from scipy.misc import comb
#    c = comb(N, 2, exact=True, repetition=False)
#    #import matplotlib.pyplot as plt
#    #plt.bar(range(0, len(p)), p, align='center')
#    if c > 2:
#        expec_real = [ 1 if round(x*c) >= 1 else 0 for x in p ]
#        EDthreshold = min(EDthreshold, expec_real.index(1) - 1)
#    return EDthreshold

def checkBCcouples_ED(any_df, inside_ShS=True, max_ED=12):
    ### USAGE:
    # Designed to be used AT IS level (i.e. inside a sample, on a single IS)
    
    ### NOTE ABOUT KWARGS:
    # inside_ShS: whether to compute EditDistances only inside
    #             ShS compartments (True) or all-vs-all (False)
    # max_ED: only data with ED <= max_ED are
    #                  stored and returned
    
    ### Fixed parameters ###########################################
    sort_bySC='descending'                                         #
    sort_byED='ascending'                                          #
    # Note: Any changes invalidate logics of                       #
    #       filterBy_randomBC_EditDistance!!!                      #
    # However, 'ascending', 'descending' or False/None are allowed #
    # if you are sure about what you're doing.                     #
    ################################################################
    
    def build_BCcouples_ED_DF(any_df, shearsite=None, max_ED=max_ED, sort_bySC=sort_bySC):
        import editdistance  # e.g.: editdistance.eval('banana', 'bahama') >>> 2L
        # Prepare Data
        DF = any_df
        # Sort data By SC if set
        if sort_bySC:
            if sort_bySC == 'ascending':
                DF = DF.sort(columns='seq_count', ascending=True)  #DF.sort(columns='seq_count', inplace=True, ascending=True)
            elif sort_bySC == 'descending':
                DF = DF.sort(columns='seq_count', ascending=False)  #DF.sort(columns='seq_count', inplace=True, ascending=False)
            else:
                print "\n[ERROR] checkBCcouples_ED wrong input! sort_bySC={sort_bySC} is not supported. Please choose 'ascending' or 'descending'.".format(sort_bySC=str(sort_bySC))
                sys.exit("\n[QUIT]\n")
        # Build BCcouples_ED_DF for randomBC couples whose edit distance is <= max_ED
        BC_list = DF['randomBC'].values
        SC_list = DF['seq_count'].values
        l = []
        append_to_l = l.append
        r = 0
        for index, row in DF.iterrows():
            r += 1
            source_BC = row['randomBC']
            source_SC = row['seq_count']
            for target_BC, target_SC in zip(BC_list[r:], SC_list[r:]):
                editDistance = int(editdistance.eval(source_BC, target_BC))
                if editDistance <= max_ED:
                    if shearsite:
                        append_to_l((source_BC, target_BC, source_SC, target_SC, editDistance, shearsite))
                    else:
                        append_to_l((source_BC, target_BC, source_SC, target_SC, editDistance))
                else:
                    continue
        # Create BCcouples_ED_DF
        col_names = None
        if shearsite:
            col_names = ['source_BC', 'target_BC', 'source_SC', 'target_SC', 'editDistance', 'shearsite']
        else:
            col_names = ['source_BC', 'target_BC', 'source_SC', 'target_SC', 'editDistance']
        BCcouples_ED_DF = pd.DataFrame.from_records(l, columns=col_names)
        return BCcouples_ED_DF
    
    # try to take required columns -> new DF
    required_columns = []
    if inside_ShS:
        required_columns = ['shearsite', 'randomBC', 'seq_count']
    else:
        required_columns = ['randomBC', 'seq_count']
    l = []
    try:
        [l.append(any_df.loc[:,c].to_frame()) for c in required_columns]
    except:
        print "\n[ERROR] checkBCcouples_ED wrong input! Columns {required_columns} are required. Given dataframe has {columns_found}.".format(required_columns=str(required_columns), columns_found=str(list(any_df)))
        sys.exit("\n[QUIT]\n")
    DF = pd.concat(l, axis=1, join='inner')
    # Shuffle DF to guarantee randomness in case of ties
    from numpy.random import permutation as shuffle_rows
    DF = DF.iloc[shuffle_rows(len(DF))]
    DF.reset_index(inplace=True, drop=True)
    
    # Create BCcouples_ED_DF from new DF
    BCcouples_ED_DF = None
    if inside_ShS:
        grouped = DF.groupby('shearsite')
        BCcouples_ED_DF = grouped.apply(lambda f: build_BCcouples_ED_DF(f, shearsite=f['shearsite'].iloc[0]))
    else:
        BCcouples_ED_DF = build_BCcouples_ED_DF(DF)
    # Sort settings for BCcouples_ED_DF
    sort_list = []
    ascending_list = []
    # by ED
    if sort_byED:
        sort_list.append('editDistance')
        if sort_byED == 'ascending':
            ascending_list.append(True)
        elif sort_byED == 'descending':
            ascending_list.append(False)
        else:
            print "\n[ERROR] checkBCcouples_ED wrong input! sort_byED={sort_byED} is not supported. Please choose 'ascending' or 'descending'.".format(sort_byED=str(sort_byED))
            sys.exit("\n[QUIT]\n")
    # by ShS if possible
    if inside_ShS:
        sort_list.append('shearsite')
        ascending_list.append(True)
    # by SC
    if sort_bySC:
        sort_list.append('source_SC')
        if sort_bySC == 'ascending':
            ascending_list.append(True)
        elif sort_bySC == 'descending':
            ascending_list.append(False)
        else:
            print "\n[ERROR] checkBCcouples_ED wrong input! sort_bySC={sort_bySC} is not supported. Please choose 'ascending' or 'descending'.".format(sort_bySC=str(sort_bySC))
            sys.exit("\n[QUIT]\n")
    # Sort BCcouples_ED_DF and reset index
    if len(sort_list) > 0:
        BCcouples_ED_DF.sort(columns=sort_list, ascending=ascending_list, inplace=True)
    BCcouples_ED_DF.reset_index(inplace=True, drop=True)
    # May be useful ###########################################################
    # TS = BCcouples_ED_DF.loc[:,'editDistance'].value_counts()               #
    # nested_array = BCcouples_ED_DF.loc[:,['shearsite', 'target_BC']].values #
    ###########################################################################
    return BCcouples_ED_DF

def filterBy_randomBC_EditDistance(any_df, inside_ShS=True, ED_rule=3):
    ### USAGE:
    # Designed to be used on a complete df (headers not required, shearsite depending on inside_ShS kwarg)
    
    ### NOTE ABOUT KWARGS:
    # inside_ShS: whether to work with EditDistances only inside
    #             ShS compartments (True) or all-vs-all (False)
    # ED_rule: can be an int or a simple callable, function of N (distinct BC);
    #          e.g." ED_rule=findEDthreshold" lead to "EDthreshold=ED_rule(len(frame))".
    #          If an int, controls will be performed. If callable, you have to be sure that an int will
    #          be returned and no controls will be performed.
    
    ### GENERAL NOTE ABOUT FUNC BEHAVIOUR:
    # exploiting checkBCcouples_ED within each group (see grouping_rule), entries with editDistance =< ED_rule
    # will be removed from the returned dataframe
    
    verbosePrint("    > Running filterBy_randomBC_EditDistance ...")
    # set grouping_rule and declare EDthreshold as None
    grouping_rule = ['association_ID', 'genomic_coordinates']
    EDthreshold = None
    # check any_df and kwargs
    verbosePrint("      * Checking data and settings ...")
    required_columns = []
    
    ######################################################################################################################################################
    
    ### >>>>>>> CODE 1 and CODE 2 ARE MUTUALLY EXCLUSIVE <<<<<<< ###
    
    ### ACTUAL VERSION ###
    ## How it works:
    ## Grouping is performed only by ['association_ID', 'genomic_coordinates']; further ShS groupby is performed performed by funcions
    ## called in the main loop. Doing this way, this function is fast and optimized for static ED_rule, e.g. D_rule = 3.
    ## However, in case of ED_rule = findEDthreshold(N), the EXPERIMENTAL VERSION should be preferred.
    ### CODE 1 #####################################################################
    if inside_ShS:                                                                 #
        required_columns = grouping_rule + ['shearsite', 'randomBC', 'seq_count']  #
    else:                                                                          #
        required_columns = grouping_rule + ['randomBC', 'seq_count']               #
    ################################################################################
    
    ### EXPERIMENTAL VERSION ###
    ## How it works:
    # Grouping by ShS right here, outside the main loop, allows better inference of ED_rule through ED_rule = findEDthreshold(N) (specific for each ShS!)
    # However, this way some operations turns out to be redundant (e.g. further ShS groupby performed by funcions called in the main loop)
    # and other coding choices may be not optimized, especially in case of static ED_rule, e.g. D_rule = 3.
    ### CODE 2 ######################################################
    #if inside_ShS:                                                 #
    #    grouping_rule = grouping_rule + ['shearsite']              #
    #required_columns = grouping_rule + ['randomBC', 'seq_count']   #
    #################################################################
    
    ######################################################################################################################################################
    
    if not set(list(any_df)).intersection(set(required_columns)) == set(required_columns):
        verbosePrint("      [Warning] Columns {required_columns} are required for running filterBy_randomBC_EditDistance (inside_ShS={inside_ShS}). Given dataframe has {columns_found}.".format(required_columns=str(required_columns), columns_found=str(list(any_df)), inside_ShS=str(inside_ShS)))
        verbosePrint("      ...skip!")
        return any_df
    if len(any_df) < 1:
        print "\n[ERROR] filterBy_randomBC_EditDistance wrong input! Given dataframe is empty!"
        sys.exit("\n[QUIT]\n")
    if type(ED_rule) is int:
        if (ED_rule >= 0) and (ED_rule < len(any_df.loc[0, 'randomBC'])):
            EDthreshold = ED_rule
        else:
            print "\n[ERROR] filterBy_randomBC_EditDistance wrong input! Must hold: 0 =< ED_rule < {n}. Your ED_rule: {ED_rule}.".format(ED_rule=str(ED_rule), n=str(len(any_df.loc[0, 'randomBC'])))
            sys.exit("\n[QUIT]\n")
    elif not callable(ED_rule):
        print "\n[ERROR] filterBy_randomBC_EditDistance wrong input! ED_rule must be a valid int or a function/callable. Your ED_rule: {ED_rule}, type={t}.".format(ED_rule=str(ED_rule), t=str(type(ED_rule)))
        sys.exit("\n[QUIT]\n")
    # Here ED_rule is an int or a callable, then EDthreshold is an int or None
    verbosePrint("      OK.")
    verbosePrint("      {n} total entries to analyze.".format(n=str(len(any_df))))
    verbosePrint("      inside_ShS: {inside_ShS}.".format(inside_ShS=str(inside_ShS)))
    if EDthreshold is not None:
        verbosePrint("      Static ED_rule: {ED_rule}. Filer out if less or equal!".format(ED_rule=str(ED_rule)))
    else:
        verbosePrint("      Dynamic ED_rule: inferred for each group.")
    # prepare data
    verbosePrint("      * Grouping data by {grouping_rule} ...".format(grouping_rule=str(grouping_rule)))
    grouped = any_df.groupby(grouping_rule)
    # Note: cannot access to frames in grouped as grouped[key] BUT you can through grouped.get_group(key)
    verbosePrint("        Done. {m} groups.".format(m=str(len(grouped))))
    # loop and filter
    filtered_df_list = []
    append_to_filtered_df_list = filtered_df_list.append
    i = 1
    # loop for inside_ShS=True case
    if inside_ShS:
        verbosePrint("      * Looping over grouped DF, within shearsite compartments, and cleaning ...")
        for key, frame in grouped:
            # print "key:", key
            # print "frame:", "\n", frame
            verbosePrint("        > {i} of {m} - processing {key} ...".format(key=str(key), i=str(i), m=str(len(grouped))))
            i += 1
            if type(ED_rule) is not int:
                EDthreshold = ED_rule(len(frame))
            if len(frame) != 1:
                verbosePrint("          {n} items, EDthreshold={EDthreshold}.".format(n=str(len(frame)), EDthreshold=str(EDthreshold)))
                BCcouples_ED_DF = checkBCcouples_ED(frame, inside_ShS=True, max_ED=EDthreshold)
                # Old code when type(shearsite) was 'str'##################################################################
                #to_remove = set([x[0]+x[1] for x in BCcouples_ED_DF[['shearsite', 'target_BC']].values])
                #filtered_frame = frame[~frame[['shearsite', 'randomBC']].apply(lambda x: x[0]+x[1] in to_remove, axis=1)]
                ###########################################################################################################
                to_remove = set([(x[0], x[1]) for x in BCcouples_ED_DF[['shearsite', 'target_BC']].values])
                filtered_frame = frame[~frame[['shearsite', 'randomBC']].apply(lambda x: (x[0], x[1]) in to_remove, axis=1)]
                append_to_filtered_df_list(filtered_frame)
            else:
                verbosePrint("          {n} item: skip!".format(n=str(len(frame))))
                append_to_filtered_df_list(frame)
    # loop for inside_ShS=False case
    else:
        verbosePrint("      * Looping over grouped DF and cleaning ...")
        for key, frame in grouped:
            # print "key:", key
            # print "frame:", "\n", frame
            verbosePrint("        > {i} of {m} - processing {key} ...".format(key=str(key), i=str(i), m=str(len(grouped))))
            i += 1
            if type(ED_rule) is not int:
                EDthreshold = ED_rule(len(frame))
            verbosePrint("          {n} items, EDthreshold={EDthreshold}.".format(n=str(len(frame)), EDthreshold=str(EDthreshold)))
            if len(frame) != 1:
                BCcouples_ED_DF = checkBCcouples_ED(frame, inside_ShS=False, max_ED=EDthreshold)
                to_remove = set(BCcouples_ED_DF['target_BC'].values)
                filtered_frame = frame[~frame['randomBC'].isin(to_remove)]
                append_to_filtered_df_list(filtered_frame)
            else:
                append_to_filtered_df_list(frame)
    # build filtered_any_df to return
    filtered_any_df = pd.concat(filtered_df_list)
    verbosePrint("        Done. {n} entries kept ({p}% discarded)".format(n=str(len(filtered_any_df)), p=str((float(len(any_df)) - float(len(filtered_any_df))) / float(len(any_df)) * 100)[:5]))
    verbosePrint("    > Done!")
    return filtered_any_df





#++++++++++++++++++++++++++++++++++++++ MAIN and TEST +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

if __name__ == "__main__":
    
    #### EXAMPLE OF USAGE (Actual version)
    filtered_df = filterBy_randomBC_EditDistance(any_df, inside_ShS=True, ED_rule=3)
    
    #### EXAMPLE OF USAGE for EXPERIMENTAL VERSION - Code 2 #####################################
    #filtered_df = filterBy_randomBC_EditDistance(df, inside_ShS=True, ED_rule=findEDthreshold) #
    #############################################################################################
    