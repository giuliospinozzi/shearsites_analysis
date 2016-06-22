# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 15:16:12 2016

@author: stefano
"""


#++++++++++++++ Requested Package(s) Import +++++++++++++++#
import os, sys
import pandas as pd
import matrix_RandomBC_globModule

#++++++++++++++++++++++ Global Vars +++++++++++++++++++++++#
#verbose = matrix_RandomBC_globModule.verbose


#++++++++++++++++++++++ Global Funcs ++++++++++++++++++++++#
verbosePrint = matrix_RandomBC_globModule.verbosePrint
humanSorted = matrix_RandomBC_globModule.humanSorted


#+++++++++++++++++++++++++++++++++++++++ FUNCTIONS ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#





def loadGzFile (path, eol='\n'):
    """
    Purpose: LOAD A GZ FILE as line list
    IN:
    path, eol
    OUT:
    file_as_list - list of file's rows without eols
    """
    import gzip
    file_content = None
    verbosePrint("\n> Loading discarded headers ...")
    verbosePrint("  > file path: {path}".format(path=str(path)))
    with gzip.open(path, 'r') as f:
        file_content = f.read()
    file_as_list = file_content.split(eol)
    verbosePrint("> headers loaded!")
    return file_as_list
    
def filterDF_byHeaders(any_df, headers_to_remove):
    verbosePrint("       > Running filterDF_byHeaders ...")
    # Check any_df
    if 'header_list' not in any_df.columns:
        print "\n[ERROR] filterDF_byHeaders wrong input! 'header_list' column is required. Given dataframe has {columns_found}.".format(columns_found=str(list(any_df)))
        sys.exit("\n[QUIT]\n")
    # Create headers_to_remove_set
    verbosePrint("         * creating the set of discarded headers ...")
    headers_to_remove_set = set(headers_to_remove)
    # Check headers_to_remove
    verbosePrint("         * checking the set of discarded headers ...")
    if len(headers_to_remove) != len(headers_to_remove_set):
        verbosePrint('''         [WARNING] filterDF_byHeaders has found duplicate headers!''')  # add here some debug stuff if useful
    # create df_header_set
    verbosePrint("         * creating the set of headers in DF ...")
    import itertools
    df_header_set = set(itertools.chain.from_iterable(any_df['header_list'].tolist()))
    # check df_header_set
    verbosePrint("         * checking the set of headers in DF ...")
    if len(df_header_set) != sum(any_df.seq_count):
        print "\n[ERROR] filterDF_byHeaders found inconsistency in DF! {n_distinct_headers} distinct headers found but the total seq_count is {tot_SC}!".format(n_distinct_headers=str(len(df_header_set)), tot_SC=str(sum(any_df.seq_count)))
        sys.exit("\n[QUIT]\n")
    # create actual_headers_to_remove_set
    verbosePrint("         * intersecting sets to get headers to remove ...")
    actual_headers_to_remove_set = headers_to_remove_set.intersection(df_header_set)
    # create filtered_DF
    filtered_DF = pd.DataFrame(columns=any_df.columns)
    # Loop over any_df and store UPDATED DATA in filtered_DF
    verbosePrint("         * Looping over DF and cleaning ...")
    n_rows = len(any_df)
    step = n_rows / 100.0
    counter = 0
    for i, r in any_df.iterrows():
        counter += 1
        if counter%int(step) == 0:
            perc = counter/step
            verbosePrint("           {perc}% ...".format(perc=str(int(perc))))
        i_header_set = set(r['header_list'])
        # case 1: all headers survived - assign row and no update
        if i_header_set.intersection(actual_headers_to_remove_set) == set():
            filtered_DF.loc[len(filtered_DF)] = r  # append row to filtered_DF
        # case 2: all headers discarded - do not assign row
        elif i_header_set.intersection(actual_headers_to_remove_set) == i_header_set:
            continue  # next row!
        # case 3: just some headers to remove - UPDATE ROW AND ASSIGN
        else:
            headers_to_keep = list(i_header_set.difference(actual_headers_to_remove_set))
            # Update header_list
            updated_r = r.set_value('header_list', headers_to_keep)
            # Update seq_count
            updated_r = r.set_value('seq_count', len(headers_to_keep))
            # Assign
            filtered_DF.loc[len(filtered_DF)] = updated_r
    verbosePrint("       > Done!")
    return filtered_DF





def filterDF_byRandomBCseqCount(any_df, SC_threshold=1, inside_ShS=True, allow_IS_loss=True, tweak_print=False):
    ### USAGE:
    # it works both inside and outside samples, use it according to what you want to do!
    
    ### NOTE ABOUT KWARGS:
    # SC_threshold: save entries > SC_threshold (put 0 to save the whole data)
    # tweak_print: False for tipical usage, or a 'barcode' string is expected, just to improve print when called by filterDF_perSample_byRandomBCseqCount
    
    ### NOTE ABOUT FUNC BEHAVIOUR:
    # 1) select required_columns based on kwargs
    # 2) groupby required_columns but seq_count
    # 3) aggregate by summing up seq_count
    # 4) identify row with seq_count > SC_threshold -> filtering_tuple_list
    #   4') if allow_IS_loss is False, filtering_tuple_list is corrected-back to reintroduce stuff whose removal causes whole IS deletion 
    # 5) compute filtered_any_df according to filtering_tuple_list
    if tweak_print:
        verbosePrint("       > {s} - Running filterDF_byRandomBCseqCount ...".format(s=str(tweak_print)))
    else:
        verbosePrint("       > Running filterDF_byRandomBCseqCount ...")
        verbosePrint("         SC_threshold: {SC_threshold}".format(SC_threshold=str(SC_threshold)))
        verbosePrint("         inside_ShS: {inside_ShS}".format(inside_ShS=str(inside_ShS)))
        verbosePrint("         allow_IS_loss: {allow_IS_loss}".format(allow_IS_loss=str(allow_IS_loss)))
    # Basic check for columns: 'randomBC' and 'seq_count';
    # if inside_ShS is True, check also for 'shearsite' column.
    # try to take required columns -> new DF
    required_columns = []
    if inside_ShS:
        required_columns = ['genomic_coordinates', 'shearsite', 'randomBC', 'seq_count']  #the order must be this
    else:
        required_columns = ['genomic_coordinates', 'randomBC', 'seq_count']  #the order must be this
    verbosePrint("         * Extracting required data: {required_columns} ...".format(required_columns=str(required_columns)))
    l = []
    try:
        [l.append(any_df.loc[:,c].to_frame()) for c in required_columns]
    except:
        print "\n[ERROR] filterDF_byRandomBCseqCount wrong input! Columns {required_columns} are required when kwarg inside_ShS={inside_ShS}. Given dataframe has {columns_found}.".format(required_columns=str(required_columns), columns_found=str(list(any_df)), inside_ShS=str(inside_ShS))
        sys.exit("\n[QUIT]\n")
    verbosePrint("         * Building extracted data ...")
    DF = pd.concat(l, axis=1, join='inner')
    verbosePrint("           Done. {n} total entries to analyze.".format(n=str(len(DF))))
    
    # Groupby and aggregate -> grouped_DF
    grouping_rule = required_columns[:-1]
    verbosePrint("         * Grouping data by {grouping_rule} and get global SC for each IS ...".format(grouping_rule=str(grouping_rule)))
    grouped = DF.groupby(grouping_rule)
    grouped_DF = grouped.sum()
    # Here grouped_DF has grouping_rule tuple as multiindex (tuple) and seq_count (the sum) as the only column
    
    # compute filtering_tuple_list: each tuple has as many items as grouping_rule, same order
    # filtering_tuple_list gives indications about which entries of any_df we have to keep
    verbosePrint("         * Flagging DF's entries to keep ...")
    filtering_tuple_list = None
    if allow_IS_loss:
        filtered_grouped_DF = grouped_DF[grouped_DF['seq_count'] > SC_threshold]
        filtering_tuple_list = filtered_grouped_DF.index.values.tolist()
    else:
        IS_orig = set(any_df['genomic_coordinates'])
        tmp_filtered_grouped_DF = grouped_DF[grouped_DF['seq_count'] > SC_threshold]
        tmp_filtered_grouped_DF.reset_index(inplace=False)
        IS_kept = set(tmp_filtered_grouped_DF.reset_index(inplace=False)['genomic_coordinates'])
        IS_lost = IS_orig.difference(IS_kept)
        if IS_lost == set():
            filtering_tuple_list = tmp_filtered_grouped_DF.index.values.tolist()
        else:
            filtering_tuple_list = tmp_filtered_grouped_DF.index.values.tolist()
            verbosePrint("         * Flagging entries whose removal would imply ISs deletion ...")
            #for genomic_coordinates in IS_lost:
            #    any_df_rows_to_restore = any_df[any_df['genomic_coordinates']==genomic_coordinates]
            #    filtering_tuple_list += any_df_rows_to_restore[grouping_rule].apply(tuple, axis=1).tolist()
            any_df_rows_to_restore = any_df[any_df['genomic_coordinates'].isin(IS_lost)]
            filtering_tuple_list += any_df_rows_to_restore[grouping_rule].apply(tuple, axis=1).tolist()
            
    # compute filtered_any_df according to filtering_tuple_list
    verbosePrint("         * Looping over DF and cleaning ...")
    filtered_any_df = any_df.copy()
    filtered_any_df['filtering_tuple'] = filtered_any_df[grouping_rule].apply(tuple, axis=1)
    filtered_any_df = filtered_any_df[filtered_any_df['filtering_tuple'].isin(filtering_tuple_list)]
    filtered_any_df.drop('filtering_tuple', axis=1, inplace=True)
    
    if tweak_print:
        verbosePrint("           Done. {n} entries kept ({p}% discarded)".format(n=str(len(filtered_any_df)), p=str((float(len(DF)) - float(len(filtered_any_df))) / float(len(DF)) * 100)[:5]))
    else:
        verbosePrint("       > Done. {n} entries kept ({p}% discarded)".format(n=str(len(filtered_any_df)), p=str((float(len(DF)) - float(len(filtered_any_df))) / float(len(DF)) * 100)[:5]))

    return filtered_any_df

def filterDF_perSample_byRandomBCseqCount(any_df, SC_threshold=1, inside_ShS=True, allow_IS_loss=True):
    verbosePrint("    > Running filterDF_perSample_byRandomBCseqCount ...")
    verbosePrint("      SC_threshold: {SC_threshold}".format(SC_threshold=str(SC_threshold)))
    verbosePrint("      inside_ShS: {inside_ShS}".format(inside_ShS=str(inside_ShS)))
    verbosePrint("      allow_IS_loss: {allow_IS_loss}".format(allow_IS_loss=str(allow_IS_loss)))
    sample_list = any_df['barcode'].unique().tolist()
    verbosePrint("      Samples found: {sample_list}.".format(sample_list=str(sample_list)))
    orig_any_df_len = len(any_df)
    verbosePrint("      Total entries to analyze: {orig_any_df_len}.".format(orig_any_df_len=str(orig_any_df_len)))
    verbosePrint("      Slicing DataFrame ...")
    l = [ any_df[any_df['barcode']==barcode] for barcode in sample_list ]
    verbosePrint("      Looping over slices ...")
    l_filt = [ filterDF_byRandomBCseqCount(x, SC_threshold=SC_threshold, inside_ShS=inside_ShS, allow_IS_loss=allow_IS_loss, tweak_print=str(x['barcode'].unique())) for x in l ]
    verbosePrint("      Building final data ...")
    filtered_any_df = pd.concat(l_filt)
    filtered_any_df.sort_index(inplace=True)
    verbosePrint("    > Done. {n} entries kept ({p}% discarded)".format(n=str(len(filtered_any_df)), p=str((float(orig_any_df_len) - float(len(filtered_any_df))) / float(orig_any_df_len) * 100)[:5]))
    return filtered_any_df





def findEDthreshold(tmp_value=3):
    ###############################
    # Still to be developed to be #
    # used as ED_rule func        #
    # Up to now it simply return  #
    # tmp_value kwarg             #
    ###############################
    
    ## MATERIAL
    ## Theor Data
    #l = [0.03167635, 0.1267054, 0.2322932, 0.2581036, 0.1935777, 0.1032414, 0.04014945, 0.01147127, 0.002389848, 0.0003540516, 0.00003540516, 0.000002145767, 0.000000005960]
    #l.reverse()
    #
    #import matplotlib.pyplot as plt
    #plt.bar(range(0, len(l)), l, align='center')
    
    return tmp_value

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
    # ED_rule: can be an int or a simple callable (no arg up to now). See "EDthreshold = ED_rule()"
    #          if an int, controls will be performed. If callable, you have to be sure that an int will
    #          be returned and no controls will be performed.
    
    ### GENERAL NOTE ABOUT FUNC BEHAVIOUR:
    # exploiting checkBCcouples_ED within each group (see grouping_rule), entries with editDistance =< ED_rule
    # will be removed from the returned dataframe
    
    verbosePrint("    > Running filterBy_randomBC_EditDistance ...")
    # set grouping_rule and declare EDthreshold as None
    grouping_rule = ['barcode', 'genomic_coordinates']
    EDthreshold = None
    # check any_df and kwargs
    verbosePrint("      * Checking df and settings ...")
    required_columns = []
    if inside_ShS:
        required_columns = grouping_rule + ['shearsite', 'randomBC', 'seq_count']
    else:
        required_columns = grouping_rule + ['randomBC', 'seq_count']
    if not set(list(any_df)).intersection(set(required_columns)) == set(required_columns):
        print "\n[ERROR] filterBy_randomBC_EditDistance (inside_ShS={inside_ShS}) wrong input! Columns {required_columns} are required. Given dataframe has {columns_found}.".format(required_columns=str(required_columns), columns_found=str(list(any_df)), inside_ShS=str(inside_ShS))
        sys.exit("\n[QUIT]\n")
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
                EDthreshold = ED_rule()
            verbosePrint("          EDthreshold: {EDthreshold}.".format(EDthreshold=str(EDthreshold)))
            verbosePrint("          {n} entries to analyze.".format(n=str(len(frame))))
            if len(frame) != 1:
                BCcouples_ED_DF = checkBCcouples_ED(frame, inside_ShS=True, max_ED=EDthreshold)
                to_remove = set([x[0]+x[1] for x in BCcouples_ED_DF[['shearsite', 'target_BC']].values])
                filtered_frame = frame[~frame[['shearsite', 'randomBC']].apply(lambda x: x[0]+x[1] in to_remove, axis=1)]
                append_to_filtered_df_list(filtered_frame)
            else:
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
                EDthreshold = ED_rule()
            verbosePrint("          EDthreshold: {EDthreshold}.".format(EDthreshold=str(EDthreshold)))
            verbosePrint("          {n} entries to analyze.".format(n=str(len(frame))))
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
    
    ### Load Association File Data
    asso_folder = "/home/stefano/Desktop/RandomBC_matrix_development/test_input/asso"
    asso_file_name = "asso.assayvalidation.lane2.tsv"
    asso_delimiter = '\t'
    # load
    import matrix_RandomBC_assoModule
    asso_dict = matrix_RandomBC_assoModule.loadAssoFile(asso_file_name, asso_folder, asso_delimiter)
    
    ### Load Data
    data_files_delimiter = '\t'
    data_files_name_filter = ".shearsites.ISfixed.LENGTHfixed.randomBC.qf100.tsv"
    develop_input_data_path = "/home/stefano/Desktop/RandomBC_matrix_development/test_input/data"
    # load
    import matrix_RandomBC_dataModule
    filtered_dir_content = matrix_RandomBC_dataModule.listDir(develop_input_data_path, name_filter=data_files_name_filter)
    verbosePrint("\n\n>>> Loading data ...\n")
    verbosePrint("> develop_input_data_path: {develop_input_data_path}".format(develop_input_data_path=str(develop_input_data_path)))
    verbosePrint("> exploited substring for data detection: '{data_files_name_filter}'".format(data_files_name_filter=str(data_files_name_filter)))
    verbosePrint("> n data files detected: {n_files}".format(n_files=str(len(filtered_dir_content))))
    verbosePrint("> data file list: {filtered_dir_content}".format(filtered_dir_content=str(filtered_dir_content)))
    verbosePrint("")
    POOL_IS_dict = {}
    POOL_alldata_dict = {}
    for path in filtered_dir_content:
        filename = str(os.path.basename(path))
        barcode = ".".join((filename.split("."))[:2])
        verbosePrint("> Processing {filename}, barcode={barcode} ...".format(filename=str(filename), barcode=str(barcode)))
        data_file_nested_list = matrix_RandomBC_dataModule.loadFile (path, data_files_delimiter)
        alldata_dict, IS_dict = matrix_RandomBC_dataModule.arrangeData(data_file_nested_list)
        POOL_IS_dict[barcode] = IS_dict
        POOL_alldata_dict[barcode] = alldata_dict
    verbosePrint("\n>>> Data Loaded!\n")
    
    #### Process Data
    import matrix_RandomBC_processingModule
    df = matrix_RandomBC_processingModule.buildDataFrame(POOL_IS_dict)
    exhaustive_df = matrix_RandomBC_processingModule.buildExhaustiveDataFrame(POOL_alldata_dict)
    
    # filtered_df = filterDF_byRandomBCseqCount(df, SC_threshold=1, inside_ShS=True, allow_IS_loss=True)
    # filtered_df = filterBy_randomBC_EditDistance(df, inside_ShS=True, ED_rule=3)
    