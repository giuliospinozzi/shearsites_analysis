# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 10:38:43 2016

@author: stefano
"""

#++++++++++++++ Requested Package(s) Import ++++++++++++++++++++++++++++++#
import sys, os
import matrix_RandomBC_globModule
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import editdistance  # e.g.: editdistance.eval('banana', 'bahama') >>> 2L


#++++++++++++++++++++++ Global Vars ++++++++++++++++++++++++++++++++++++++#
verbose = matrix_RandomBC_globModule.verbose


#++++++++++++++++++++++ Global Funcs +++++++++++++++++++++++++++++++++++++#
verbosePrint = matrix_RandomBC_globModule.verbosePrint
humanSorted = matrix_RandomBC_globModule.humanSorted


#+++++++++++++++++++++++++++++++++++++++ FUNCTIONS +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#



def checkNucleotideBalancing(any_df):
    # try to take required columns -> new DF
    required_columns = ['randomBC', 'seq_count']
    l = []
    try:
        [l.append(any_df.loc[:,c].to_frame()) for c in required_columns]
    except:
        print "\n[ERROR] checkNucleotideBalancing wrong input! Columns {required_columns} are required. Given dataframe has {columns_found}.".format(required_columns=str(required_columns), columns_found=str(list(any_df)))
        sys.exit("\n[QUIT]\n")
    DF = pd.concat(l, axis=1, join='inner')
    # operate on DF to create BC_DF: rows are RandomBC exploded by seq_count, splitted in columns on each nucleotide
    BC_DF_rowlist = []
    BC_DF_rowlist_append = BC_DF_rowlist.append
    for index, columns in DF.iterrows():
        [BC_DF_rowlist_append(list(columns['randomBC'])) for times in range(columns['seq_count'])]
    BC_DF =pd.DataFrame(BC_DF_rowlist)
    # get value counts for each columns of BC_DF
    l = []
    for n in list(BC_DF):  # columns are nucleotides positions 0-indexed
        l.append(BC_DF.loc[:,n].value_counts())
    # create nucleotidesCount_DF: nucleotide types are row indexes,
    # nucleotide positions are column indexes, values are counts.
    nucleotidesCount_DF = pd.concat(l, axis=1)
    #######################################################################################
    # eg of output usage: PLOT ----> nucleotidesCount_DF.T.plot(kind='bar', stacked=True) #
    # for better graphics please call plotNucleotideBalancing(nucleotidesCount_DF)        #
    #######################################################################################
    return nucleotidesCount_DF
    
def plotNucleotideBalancing(nucleotidesCount_DF, title='PILED-UP RANDOM-BARCODES', stacked_bar=True, show_live=True, export=""):
    # Note: nucleotidesCount_DF from checkNucleotideBalancing
    # Note: export is supposed to be a complete-and-valid path.
    #       False, None, empty-string means "no export"
    
    # Rename columns of nucleotidesCount_DF
    columns_relabelling = {}
    [columns_relabelling.update({i:i+1}) for i in range(len(list(nucleotidesCount_DF)))]
    nucleotidesCount_DF = nucleotidesCount_DF.rename(columns=columns_relabelling)
    # Set up interactive mode (plot pop-upping)
    if show_live:
        plt.ion()
    else:
        plt.ioff()
    # Prepare plot environment
    fig = plt.figure() # Create matplotlib figure
    plt.title(title)
    ax = fig.add_subplot(111) # Create matplotlib axes
    ax2 = ax.twinx() # Create another axes that shares the same x-axis as ax
    # Prepare axis
    ax.set_xlabel('base position')
    ax.set_ylabel('N of nucleotides per base')
    ax.set_ylim([0, max(nucleotidesCount_DF.sum())])
    ax2.set_ylabel('% of nucleotides per base')
    ax2.set_ylim([0,100])
    # Plot data
    nucleotidesCount_DF.T.plot(kind='bar', stacked=stacked_bar, ax=ax, position=1)
    if show_live:
        plt.show()
    # Save plot
    if export:
        export = os.path.normpath(export)
        plt.savefig(export)
    # close
    if not show_live:
        plt.close()



def checkRandomBCfrequency(any_df):
    # try to take required columns -> new DF
    required_columns = ['randomBC', 'seq_count']
    l = []
    try:
        [l.append(any_df.loc[:,c].to_frame()) for c in required_columns]
    except:
        print "\n[ERROR] checkRandomBCfrequency wrong input! Columns {required_columns} are required. Given dataframe has {columns_found}.".format(required_columns=str(required_columns), columns_found=str(list(any_df)))
        sys.exit("\n[QUIT]\n")
    DF = pd.concat(l, axis=1, join='inner')
    # operate on DF to create distinctBC_DF
    # create distinctBC_DF: distinct randomBC are row indexes,
    # seq_count is the only column, values are randomBC counts.
    distinctBC_DF = pd.pivot_table(DF, columns='randomBC', values='seq_count', aggfunc=sum).to_frame()
    #######################################################################################
    # eg of output usage: PLOT ----> distinctBC_DF.plot(kind='bar')                       #
    # however this is not usually feasible due to the large amount of distinct BC         #
    # for typical usage please call plotRandomBCfrequency(distinctBC_DF)                  #
    #######################################################################################
    return distinctBC_DF

def plotRandomBCfrequency(distinctBC_DF, title='RANDOM-BARCODE FREQUENCY', show_top_ranked=10, show_live=True, export=""):
    # Note: distinctBC_DF from checkRandomBCfrequency
    # Note: export is supposed to be a complete-and-valid path.
    #       False, None, empty-string mean "no export"
    # Note: show_top_ranked is supposed to be an int.
    #       False, None, '0' mean "do not print top ranked barcodes"
    
    # Sort data
    distinctBC_DF = pd.DataFrame.sort(distinctBC_DF, columns='seq_count', ascending=False)
    # Prepare text data if show_top_ranked
    top_represented_rBC = None
    if show_top_ranked:
        N = show_top_ranked
        top_represented_rBC = "TOP-{N} RANDOM-BARCODES:".format(N=str(N))
        for i in range(N):
            top_represented_rBC += "\n {i}) {rBC}".format(i=str(i+1), rBC=str(list(distinctBC_DF.index.values)[i]))
    # Prepare data to plot
    distinctBC_DF = distinctBC_DF.reset_index()  # distinctBC_DF is sorted so reset_index() yields the ranking as index!
    distinctBC_DF = distinctBC_DF.loc[:,'seq_count']  # take only data to plot
    # Set up interactive mode (plot pop-upping)
    if show_live:
        plt.ion()
    else:
        plt.ioff()
    # Prepare plot environment
    fig = plt.figure() # Create matplotlib figure
    plt.title(title)
    ax = fig.add_subplot(111) # Create matplotlib axes
    ax2 = ax.twinx() # Create another axes that shares the same x-axis as ax
    # Prepare axis
    y2_lim = 100 * float(distinctBC_DF.max()) / float(distinctBC_DF.sum())
    ax.set_xlabel('distinct barcodes ranked by count (N={N})'.format(N=str(len(list(distinctBC_DF)))))
    ax.set_ylabel('count (max={M}, min={m}, #classes={n})'.format(M=str(distinctBC_DF.max()), m=str(distinctBC_DF.min()), n=str(len(distinctBC_DF.unique()))))
    ax.set_ylim([0, distinctBC_DF.max()])
    ax2.set_ylabel('%')
    ax2.set_ylim([0,y2_lim])
    # Plot text data if show_top_ranked
    if show_top_ranked:
        ax.text(1.0/2.0*len(list(distinctBC_DF)), 1.0/2.0*distinctBC_DF.max(), top_represented_rBC, bbox={'facecolor':'blue', 'alpha':0.3, 'pad':10})
    # Plot data
    distinctBC_DF.plot(kind='area', ax=ax)
    if show_live:
        plt.show()
    # Save plot
    if export:
        export = os.path.normpath(export)
        plt.savefig(export)
    # close
    if not show_live:
        plt.close()



def checkEditDistance(any_df):
    # Note: 'randomBC' column required. 'shearsite' column as user whishes. 'genomic_coordinates' column exploited just for warning,
    #       if present. Any other column in any_df is ignored and cannot provide any further splitting/hierarchy on data.
    # Note: This function is DESIGNED TO BE USED ON DATA FROM A SINGLE INTEGRATION SITE, ALONG ITS SHEARSITES;
    #       This is because within an IS, for each shersite, randomBCs are supposed to be UNIQUE (they usually have a seq_count!).
    #       HOWEVER, THIS FUNCTION ALWAYS WORKS:
    #          - if more than one IS is present, data will be still processed as belonging to one unique IS
    #          - if no shearsites are present, data will be processed as belonging to one unique shear site of nominal length '0'
    #          - when any (or all) among above cases happen, uniqueness of randomBCs will be forced, due to groupby('shearsite')['randomBC'].apply(lambda x: np.unique(x.tolist()))
    
    def editDistanceMatrix (randomBC_list, shearsite_label=""):
        # Note: randomBC_list must be a list of unique BCs, supposed to be unique
        # Note: shearsite_label add a prefix to columns/rows labels/indexes; 
        #       with shearsite_label="" they remain just the BC's in randomBC_list
        sorted_unique_randomBC_list = sorted(list(set(randomBC_list)))
        if sorted_unique_randomBC_list != sorted(list(randomBC_list)):
            print "\n[ERROR] checkEditDistance (through editDistanceMatrix) found duplicate randomBCs!"
            sys.exit("\n[QUIT]\n")
        l = []
        l_append = l.append
        shearsite_label = str(shearsite_label)
        for col_rBC in sorted_unique_randomBC_list:
            key_col_rBC = shearsite_label+col_rBC
            d = {key_col_rBC: []}
            d_append = d[key_col_rBC].append
            for row_rBC in sorted_unique_randomBC_list:
                d_append(editdistance.eval(col_rBC, row_rBC))
            key_row_rBC = [shearsite_label+row_rBC for row_rBC in sorted_unique_randomBC_list]
            l_append(pd.DataFrame(data=d, index=key_row_rBC))
        ED_matrix = pd.concat(l, axis=1, join='inner')
        return ED_matrix
    
    # data collector
    l = []
    # try to take required columns -> l
    required_columns = ['randomBC']
    try:
        [l.append(any_df.loc[:,c].to_frame()) for c in required_columns]
    except:
        print "\n[ERROR] checkEditDistance wrong input! Columns {required_columns} are required. Given dataframe has {columns_found}.".format(required_columns=str(required_columns), columns_found=str(list(any_df)))
        sys.exit("\n[QUIT]\n")
    
    # check for correct usage
    try:
        IS_set = set(any_df.loc[:,'genomic_coordinates'].tolist())
        if len(IS_set) != 1:
            print "\n[WARNING] checkEditDistance found data from more than one IS found: {IS_list}. However, all the data will be processed as belonging to one unique IS!".format(IS_list=str(humanSorted(list(IS_set))))
    except:
        print "\n[WARNING] checkEditDistance cannot perform IS control on input data ('genomic_coordinates' not found). Thus, all the data will be processed as belonging to one unique IS!"
    
    # try to take other required columns, and fix if absent -> l
    try:
        l.append(any_df.loc[:,'shearsite'].to_frame())
    except:
        print "\n[WARNING] checkEditDistance can't find 'shearsite' data. All the data will be processed as belonging to one unique shear site of nominal length '0'."
        l.append(pd.DataFrame(data={'shearsite':['0']*len(any_df.loc[:,'randomBC'])}))
    
    # Prepare DF: concat(l) --> DF
    DF = pd.concat(l, axis=1, join='inner')
    # Group (append!) randomBCs by shearsite on DF
    DF = DF.groupby('shearsite')['randomBC'].apply(lambda x: np.unique(x.tolist()))
    # Here DF is returned as Series with Name=randomBC, Index=shearsite and values are lists of unique randomBC (nb: within the same list)
    # (np.unique is redundant in the best case use. Otherwise is forced. See note above)
    
    # Collect edit distance matrixe dataframes in l
    l = []
    l_append =l.append
    for shearsite_label, randomBC_list in DF.iteritems():
        # here kwarg shearsite_label is fundamental to l_append dataframes with mutually disjoint columns/rows labels/indexes
        # and then exploit 'join' as a quick way to reshape data as a whole
        l_append(editDistanceMatrix(randomBC_list, shearsite_label=shearsite_label+"_", ))
    # Reshape data in a unique dataframe using multiple joins, exploiting unique columns and rows labels and indexes!
    editDistance_DF = pd.DataFrame().join(l, how='outer')
    return editDistance_DF




#++++++++++++++++++++++++++++++++++++++ MAIN and TEST ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

if __name__ == "__main__":
    
    ### Load Association File Data
    asso_folder = "/home/stefano/Desktop/RandomBC_matrix_development/test_input/asso"
    asso_file_name = "asso.assayvalidation.lane1.tsv"
    asso_delimiter = '\t'
    # load
    import matrix_RandomBC_assoModule
    asso_dict = matrix_RandomBC_assoModule.loadAssoFile(asso_file_name, asso_folder, asso_delimiter)
    
    ### Load Data
    data_files_delimiter = '\t'
    data_files_name_filter = ".randomBC.tsv"
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
    # Up to now exhaustive_df is just like df with 'header' column added (header list).
    # Keep in mind that this structure is ongoing and more columns will may appear. Let them be implicitly supported!!
    # Check advancements in matrix_RandomBC_processingModule.buildExhaustiveDataFrame
    
    #### Test Functions in this module
    #overall_nucleotidesCount_DF = checkNucleotideBalancing(exhaustive_df)  # or = checkNucleotideBalancing(df)
    #plotNucleotideBalancing(overall_nucleotidesCount_DF, title='[DEBUG] PILED-UP SEQUENCES', stacked_bar=True, show_live=True, export=os.path.join(os.getcwd(), "test_output", "debug_checkNucleotidesBalancing.pdf"))
    #overall_distinctBC_DF = checkRandomBCfrequency(exhaustive_df)
    #plotRandomBCfrequency(overall_distinctBC_DF, title='[DEBUG] RANDOM-BARCODE FREQUENCY', show_top_ranked=10, show_live=True, export=os.path.join(os.getcwd(), "test_output", "debug_checkRandomBCfrequency.pdf"))
    # Any df is accepted so the contents of each plot is the whole input DF.
    # title and export kwargs allows you to put the proper label to your contents!