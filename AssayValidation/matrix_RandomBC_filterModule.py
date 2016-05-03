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
#humanSorted = matrix_RandomBC_globModule.humanSorted


#+++++++++++++++++++++++++++++++++++++++ FUNCTIONS ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

def buildInputPath(ground_dir, DISEASE, PATIENT):
    '''
    Purpose:
    Build and return the standard input path (thought as containing all the input data for filtering operations),
    starting from ground_dir
    IN:
    ground_dir (abs valid path), DISEASE, PATIENT (strings)
    OUT: INDIR = $ground_dir/$DISEASE/$PATIENT/'quality' (checked)
    '''
    # Check ground_dir
    try:
        ground_dir = os.path.normpath(ground_dir)
    except Exception, err_message:
        print "\n[ERROR] ground_dir='{ground_dir}' is not formatted as valid path!".format(ground_dir=str(ground_dir))
        print "os.path.normpath returned: ", err_message
        sys.exit("\n[QUIT]\n")
    if os.path.isfile(ground_dir):
        print "\n[ERROR] ground_dir must be a folder path, not a file! Your input: ground_dir='{ground_dir}'".format(ground_dir=str(ground_dir))
        sys.exit("\n[QUIT]\n")
    if not os.path.exists(ground_dir):
        print "\n[ERROR] ground_dir='{ground_dir}' does not exist!".format(ground_dir=str(ground_dir))
        sys.exit("\n[QUIT]\n")
    if not os.access(ground_dir, os.R_OK):
        print "\n[ERROR] You have not read permission in ground_dir='{ground_dir}'".format(ground_dir=str(ground_dir))
        sys.exit("\n[QUIT]\n")
    # Build INDIR
    BASEDIR = os.path.normpath(os.path.join(ground_dir, DISEASE, PATIENT))
    INDIR = os.path.normpath(os.path.join(BASEDIR, "quality"))
    # Check INDIR
    if not os.path.exists(INDIR):
        print "\n[ERROR] INDIR='{INDIR}' just build does not exist!".format(INDIR=str(INDIR))
        sys.exit("\n[QUIT]\n")
    if not os.access(INDIR, os.R_OK):
        print "\n[ERROR] You have not read permission in INDIR='{INDIR}'".format(INDIR=str(INDIR))
        sys.exit("\n[QUIT]\n")
    return INDIR

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
    verbosePrint("\n> Running filterDF_byHeaders ...")
    # Check any_df
    if 'header_list' not in any_df.columns:
        print "\n[ERROR] filterDF_byHeaders wrong input! 'header_list' column is required. Given dataframe has {columns_found}.".format(columns_found=str(list(any_df)))
        sys.exit("\n[QUIT]\n")
    # Create headers_to_remove_set
    verbosePrint("> creating the set of discarded headers ...")
    headers_to_remove_set = set(headers_to_remove)
    # Check headers_to_remove
    verbosePrint("> checking the set of discarded headers ...")
    if len(headers_to_remove) != len(headers_to_remove_set):
        verbosePrint('''[WARNING] filterDF_byHeaders has found duplicate headers!''')  # add here some debug stuff if useful
    # create df_header_set
    verbosePrint("> creating the set of headers in DF ...")
    import itertools
    df_header_set = set(itertools.chain.from_iterable(any_df['header_list'].tolist()))
    # check df_header_set
    verbosePrint("> checking the set of headers in DF ...")
    if len(df_header_set) != sum(any_df.seq_count):
        print "\n[ERROR] filterDF_byHeaders found inconsistency in DF! {n_distinct_headers} distinct headers found but the total seq_count is {tot_SC}!".format(n_distinct_headers=str(len(df_header_set)), tot_SC=str(sum(any_df.seq_count)))
        sys.exit("\n[QUIT]\n")
    # create actual_headers_to_remove_set
    verbosePrint("> intersecting sets to get headers to remove ...")
    actual_headers_to_remove_set = headers_to_remove_set.intersection(df_header_set)
    # create filtered_DF
    filtered_DF = pd.DataFrame(columns=any_df.columns)
    # Loop over any_df and store UPDATED DATA in filtered_DF
    verbosePrint("> Looping over DF and cleaning ...")
    n_rows = len(any_df)
    step = n_rows / 100.0
    counter = 0
    for i, r in any_df.iterrows():
        counter += 1
        if counter%int(step) == 0:
            perc = counter/step
            verbosePrint("  > {perc}% ...".format(perc=str(int(perc))))
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
    verbosePrint("> filtered_DF built!")
    return filtered_DF



def filterDF_byRandomBCseqCount(any_df, SC_threshold=1, inside_ShS=True, allow_IS_loss=True):
    ### WARNING:
    # it works in/out-side samples, use it according to what you want to do!
    
    ### NOTE ABOUT KWARGS:
    # SC_threshold: save entries > SC_threshold (put 0 to save the whole data)
    
    ### NOTE ABOUT FUNC BEHAVIOUR:
    # 1) select required_columns based on kwargs
    # 2) groupby required_columns but seq_count
    # 3) aggregate by summing up seq_count
    # 4) identify row with seq_count > SC_threshold -> filtering_tuple_list
    #   4') if allow_IS_loss is False, filtering_tuple_list is corrected-back to reintroduce stuff whose removal causes whole IS deletion 
    # 5) compute filtered_any_df according to filtering_tuple_list
    
    
    # Basic check for columns: 'randomBC' and 'seq_count';
    # if inside_ShS is True, check also for 'shearsite' column.
    # try to take required columns -> new DF
    required_columns = []
    if inside_ShS:
        required_columns = ['genomic_coordinates', 'shearsite', 'randomBC', 'seq_count']  #the order must be this
    else:
        required_columns = ['genomic_coordinates', 'randomBC', 'seq_count']  #the order must be this
    l = []
    try:
        [l.append(any_df.loc[:,c].to_frame()) for c in required_columns]
    except:
        print "\n[ERROR] filterDF_byRandomBCseqCount wrong input! Columns {required_columns} are required when kwarg inside_ShS={inside_ShS}. Given dataframe has {columns_found}.".format(required_columns=str(required_columns), columns_found=str(list(any_df)), inside_ShS=str(inside_ShS))
        sys.exit("\n[QUIT]\n")
    DF = pd.concat(l, axis=1, join='inner')
    
    # Groupby and aggregate -> grouped_DF
    grouping_rule = required_columns[:-1]
    grouped = DF.groupby(grouping_rule)
    grouped_DF = grouped.sum()
    # Here grouped_DF has grouping_rule tuple as multiindex (tuple) and seq_count (the sum) as the only column
    
    # compute filtering_tuple_list: each tuple has as many items as grouping_rule, same order
    # filtering_tuple_list gives indications about which entries of any_df we have to keep
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
            for genomic_coordinates in IS_lost:
                any_df_rows_to_restore = any_df[any_df['genomic_coordinates']==genomic_coordinates]
                filtering_tuple_list += any_df_rows_to_restore[grouping_rule].apply(tuple, axis=1).tolist()
                
    # compute filtered_any_df according to filtering_tuple_list
    filtered_any_df = any_df.copy()
    filtered_any_df['filtering_tuple'] = filtered_any_df[grouping_rule].apply(tuple, axis=1)
    filtered_any_df = filtered_any_df[filtered_any_df['filtering_tuple'].isin(filtering_tuple_list)]
    filtered_any_df.drop('filtering_tuple', axis=1, inplace=True)
    
    return filtered_any_df

#++++++++++++++++++++++++++++++++++++++ MAIN and TEST +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#


