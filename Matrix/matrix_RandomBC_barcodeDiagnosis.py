# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 10:38:43 2016

@author: stefano
"""

#++++++++++++++ Requested Package(s) Import ++++++++++++++++++++++++++++++#
import sys, os
import matrix_RandomBC_globModule
import pandas as pd
import matplotlib.pyplot as plt
#import editdistance  # e.g.: editdistance.eval('banana', 'bahama') >>> 2L


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
    #######################################################################################
    return nucleotidesCount_DF
    
def plotNucleotideBalancing(nucleotidesCount_DF, title='PILED-UP SEQUENCES', stacked_bar=True, show_live=True, export=""):
    # Note: nucleotidesCount_DF from checkNucleotideBalancing
    # Note: export is supposed to be a complete-and-valid path.
    #       False, None, empty-string means "no export"
    
    plt.close() # it's a try
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
    # Save plot
    if export:
        export = os.path.normpath(export)
        plt.savefig(export)


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
    return distinctBC_DF

def plotRandomBCfrequency(distinctBC_DF, show_live=True):
    distinctBC_DF = pd.DataFrame.sort(distinctBC_DF, columns='seq_count', ascending=False)
    distinctBC_DF = distinctBC_DF.reset_index()
    distinctBC_DF = distinctBC_DF.loc[:,'seq_count']
    plt.close() # it's a try
    # Set up interactive mode (plot pop-upping)
    if show_live:
        plt.ion()
    else:
        plt.ioff()
    distinctBC_DF.plot(kind='bar')
    # sistemare bene la grafica

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
    