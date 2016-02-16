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
import seaborn.apionly as sns  # this way mpl defaults are kept

from math import ceil

#++++++++++++++++++++++ Global Vars ++++++++++++++++++++++++++++++++++++++#
verbose = matrix_RandomBC_globModule.verbose


#++++++++++++++++++++++ Global Funcs +++++++++++++++++++++++++++++++++++++#
verbosePrint = matrix_RandomBC_globModule.verbosePrint
humanSorted = matrix_RandomBC_globModule.humanSorted


#+++++++++++++++++++++++++++++++++++++++ FUNCTIONS +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#



def checkNucleotideBalancing(any_df, how='by_seq_count', N_warn=False):
    # Check 'how' kwarg (it change the logics!):
    #    how='by_seq_count' explodes 'randomBC' by their 'seq_count'
    #    how='simple' can be thought as 'by_seq_count' where, for each cell of 'randomBC' column, the seq_count is imposed to be ONE
    #    how='distinct' takes all the distinct randomBCs ONCE.
    # N_warn = True produces warnings (if verbose) if randomBCs show 'N'/'n' nucleotides in their sequences
    supported_how = ['by_seq_count', 'simple', 'distinct']
    if how not in supported_how:
        print "\n[ERROR] checkNucleotideBalancing wrong input! Supported value for 'how' kwarg are: {supported_how}. Your input: {how}.".format(supported_how=str(supported_how), how=str(how))
        sys.exit("\n[QUIT]\n")
    # prepare l with selected columns to build a new DF
    l = []
    if how == 'by_seq_count':
        # try to take required columns -> append to l -> new DF
        required_columns = ['randomBC', 'seq_count']
        try:
            [l.append(any_df.loc[:,c].to_frame()) for c in required_columns]
        except:
            print "\n[ERROR] checkNucleotideBalancing wrong input! Columns {required_columns} are required. Given dataframe has {columns_found}.".format(required_columns=str(required_columns), columns_found=str(list(any_df)))
            sys.exit("\n[QUIT]\n")
    elif (how == 'simple' or  how == 'distinct'):
        # try to take required columns + fake seq_count column filled with one's -> append to l -> new DF
        required_columns = ['randomBC']
        try:
            [l.append(any_df.loc[:,c].to_frame()) for c in required_columns]
        except:
            print "\n[ERROR] checkNucleotideBalancing wrong input! Columns {required_columns} are required. Given dataframe has {columns_found}.".format(required_columns=str(required_columns), columns_found=str(list(any_df)))
            sys.exit("\n[QUIT]\n")
        # append to l a fake seq_count column filled with one's
        l.append(pd.DataFrame(data=[1]*len(any_df.loc[:,'randomBC']), index=any_df.index, columns=['seq_count']))
    # add column 'header_list' if N_warn is requested
    if N_warn:
        header_list = True
        try:
           l.append(any_df.loc[:,'header_list'].to_frame())
        except:
            header_list = False
    # Create DF concatenating df's in l
    DF = pd.concat(l, axis=1, join='inner')
    if how == 'distinct':
        DF.drop_duplicates(subset='randomBC', inplace=True)
    # operate on DF to create BC_DF: rows are RandomBC (distinct if how == 'distinct') exploded by seq_count (1's if how == 'simple' or 'distinct'), splitted in columns on each nucleotide
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
    if N_warn:
        nBC_list = []
        nBC_index_list = []
        nucleotides_found = list(nucleotidesCount_DF.index.values)
        if (('N' in nucleotides_found) or ('n' in nucleotides_found)):
            verbosePrint('''[WARNING] checkNucleotideBalancing found "N's" in sequences!''')
            for i, r in any_df.loc[:,'randomBC'].to_frame().iterrows():
                if (('N' in r['randomBC']) or ('n' in r['randomBC'])):
                    nBC_list.append(r['randomBC'])
                    nBC_index_list.append(i)
            verbosePrint('''          * randomBC with N's: {nBC_list}'''.format(nBC_list=str(nBC_list)))
            if header_list:
                header_list = []
                for i, r in any_df.loc[nBC_index_list,'header_list'].to_frame().iterrows():
                    header_list.append(r['header_list'])
                if how == 'distinct':
                    verbosePrint('''          * example headers (not exhaustive!): {header_list}'''.format(header_list=str(header_list)))
                else:
                    verbosePrint('''          * header list: {header_list}'''.format(header_list=str(header_list)))
            else:
                verbosePrint('''          * Header list is not available.''')
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
    fig = plt.figure(figsize=(13,13)) # Create matplotlib figure
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




def checkRandomBCoccurrency(any_df):
    # try to take required columns -> new DF
    required_columns = ['randomBC', 'seq_count']
    l = []
    try:
        [l.append(any_df.loc[:,c].to_frame()) for c in required_columns]
    except:
        print "\n[ERROR] checkRandomBCoccurrency wrong input! Columns {required_columns} are required. Given dataframe has {columns_found}.".format(required_columns=str(required_columns), columns_found=str(list(any_df)))
        sys.exit("\n[QUIT]\n")
    DF = pd.concat(l, axis=1, join='inner')
    # operate on DF to create distinctBC_DF
    # create distinctBC_DF: distinct randomBC are row indexes,
    # seq_count is the only column, values are randomBC counts.
    distinctBC_DF = pd.pivot_table(DF, columns='randomBC', values='seq_count', aggfunc=sum).to_frame()
    #######################################################################################
    # eg of output usage: PLOT ----> distinctBC_DF.plot(kind='bar')                       #
    # however this is not usually feasible due to the large amount of distinct BC         #
    # for typical usage please call plotRandomBCoccurrency(distinctBC_DF)                 #
    #######################################################################################
    return distinctBC_DF

def plotRandomBCoccurrency(distinctBC_DF, title='RANDOM-BARCODE OCCURRENCIES', show_top_ranked=10, annot=True, show_live=True, export=""):
    # Note: distinctBC_DF from checkRandomBCoccurrency
    # Note: export is supposed to be a complete-and-valid path.
    #       False, None, empty-string mean "no export"
    # Note: show_top_ranked is supposed to be an int.
    #       False, None, '0' mean "do not print top ranked barcodes"
    
    # Sort data
    distinctBC_DF = pd.DataFrame.sort(distinctBC_DF, columns='seq_count', ascending=False)
    # Prepare text data if show_top_ranked
    top_represented_rBC = None
    if show_top_ranked:
        N = min(show_top_ranked, len(distinctBC_DF))
        top_represented_rBC = "TOP-{N} RANDOM-BARCODES BY SEQ-COUNT:".format(N=str(N))
        for i in range(N):
            top_represented_rBC += "\n {i}) {rBC} (seq_count={sc}, {p}%)".format(i=str(i+1), rBC=str(list(distinctBC_DF.index.values)[i]), sc=str(int(distinctBC_DF['seq_count'][i])), p=str(round(list(distinctBC_DF['seq_count'])[i]*100.0/distinctBC_DF.sum(),2)))
    # Prepare data to plot as Series
    distinctBC_DF = distinctBC_DF.reset_index()  # distinctBC_DF is sorted so reset_index() yields the ranking as index!
    distinctBC_DF = distinctBC_DF.loc[:,'seq_count']  # take only data to plot: distinctBC_DF is a Series from now on
    # Set up interactive mode (plot pop-upping)
    if show_live:
        plt.ion()
    else:
        plt.ioff()
    # Prepare plot environment
    fig = plt.figure(figsize=(13,13)) # Create matplotlib figure
    plt.title(title)
    ax = fig.add_subplot(111) # Create matplotlib axes
    ax2 = ax.twinx() # Create another axes that shares the same x-axis as ax
    # Prepare axis
    y2_lim = 100 * float(distinctBC_DF.max()) / float(distinctBC_DF.sum())
    ax.set_ylabel('occurrencies (max={M}, min={m}, total={n})'.format(M=str(distinctBC_DF.max()), m=str(distinctBC_DF.min()), n=str(distinctBC_DF.sum())))
    ax.set_ylim([0, distinctBC_DF.max()])
    ax2.set_ylabel('%')
    ax2.set_ylim([0,y2_lim])
    # Plot text data if show_top_ranked
    if show_top_ranked:
        ax.text(1.0/2.0*len(list(distinctBC_DF)), 1.0/2.0*distinctBC_DF.max(), top_represented_rBC, bbox={'facecolor':'blue', 'alpha':0.3, 'pad':10})
    # Plot data
    distinctBC_DF.plot(kind='area', ax=ax)
    # Optional: annotate x axis with counts and percentages
    if annot:
        # Label the raw counts and the percentages below the x-axis
        cumulative_count = 0
        cumulative_percent = 0
        control = 0
        for x, count in distinctBC_DF.iteritems():
            cumulative_count += count
            rule = round((float(cumulative_count) / distinctBC_DF.sum())/2,1)
            if  rule != control:
                control = rule
                # Label with cumulative count
                cumulative_count_label = str(int(cumulative_count))
                ax.annotate(cumulative_count_label, xy=(x, 0), xycoords=('data', 'axes fraction'), xytext=(0, -18), textcoords='offset points', va='top', ha='center')
                # Label with cumulative percentages
                cumulative_percent = "{:.1f}".format(100 * float(cumulative_count) / distinctBC_DF.sum()) + "%"
                ax.annotate(cumulative_percent, xy=(x, 0), xycoords=('data', 'axes fraction'), xytext=(0, -32), textcoords='offset points', va='top', ha='center')
        # Refine axis
        ax.set_xlabel("distinct barcodes ranked by count (N={N})\ncumulative counts\ncumulative percentages".format(N=str(len(list(distinctBC_DF)))))
        ax.xaxis.set_label_coords(0.8, -0.03)
        # Give more room at the bottom of the plot
        plt.subplots_adjust(bottom=0.2)
    else:
        # standard X axis label
        ax.set_xlabel('distinct barcodes ranked by count (N={N})'.format(N=str(len(list(distinctBC_DF)))))
    # Show live
    if show_live:
        plt.show()
    # Save plot
    if export:
        export = os.path.normpath(export)
        plt.savefig(export)
    # close
    if not show_live:
        plt.close()




def checkEditDistance(any_df, all_combinations=False):
    # Note about any_df: 'randomBC' column required. 
    #                    'shearsite' column is optional. It splits the data in groups of equal 'shersite'.
    #                    'genomic_coordinates' is 'very' optional and may just produce warnings.
    #                    Any other column in any_df is just ignored and cannot provide any further splitting/hierarchy on data.
    # Note about all_combinations: if 'shearsite' column is in any_df, it states whether to compute edit distances only within
    #                              shearsite groups (False) and nan elsewere or produce a complete matrix (True)
    # Note: This function is DESIGNED TO BE USED ON DATA FROM A SINGLE INTEGRATION SITE, ALONG ITS SHEARSITES;
    #       This is because within an IS, for each shersite, randomBCs are supposed to be UNIQUE (they usually have a seq_count!).
    #       HOWEVER, THIS FUNCTION ALWAYS WORKS:
    #          - even if more than one IS is present, or 'genomic_coordinates' is just not available, data will be still processed as belonging to one unique IS (a warn is printed if verbose)
    #          - if no shearsites are present, data will be processed as belonging to one unique shear site of nominal length '0' (a warn is printed if verbose)
    #          - IMPORTANT: when any (or all) among above cases happen, uniqueness of randomBCs will be forced, due to groupby('shearsite')['randomBC'].apply(lambda x: np.unique(x.tolist()))
        
    def editDistanceMatrix (seqList1, seqList2, groupLabel1="", groupLabel2=""):
        # return a matrix with col from seqList1 and row from seqList2
        # groupLabel kwargs add a fixed string to row and col names
        sorted_unique_seqList1 = sorted(list(set(seqList1)))
        sorted_unique_seqList2 = sorted(list(set(seqList2)))
        l = []
        l_append = l.append
        for seq_col in sorted_unique_seqList1:
            col_key = groupLabel1+seq_col
            d = {col_key: []}
            d_append = d[col_key].append
            for seq_row in sorted_unique_seqList2:
                d_append(editdistance.eval(seq_col, seq_row))
            row_keys = [groupLabel2+seq_row for seq_row in seqList2]
            l_append(pd.DataFrame(data=d, index=row_keys))
        ED_matrix = pd.concat(l, axis=1, join='inner')
        return ED_matrix
        
    def buildGroupLabel (shearsite, prefix="", suffix="", sep="_", ignore=0):
        # ignore=0 to waste fake shearsite=0
        if str(shearsite) == str(ignore):
            if ((not prefix) and (not suffix)):
                return ""
            elif ((not prefix) or (not suffix)):
                return prefix+suffix+sep
            else:
                return sep.join([prefix,suffix,""])
        else:
            if ((not prefix) and (not suffix)):
                return str(shearsite)+sep
            elif (not prefix):
                return sep.join([str(shearsite),suffix,""])
            elif (not suffix):
                return sep.join([prefix,str(shearsite),""])
            else:
                return sep.join([prefix,str(shearsite),suffix,""])
    # data collector
    l = []
    # try to take required columns -> l
    required_columns = ['randomBC']
    try:
        [l.append(any_df.loc[:,c].to_frame()) for c in required_columns]
    except:
        print "\n[ERROR] checkEditDistance wrong input! Columns {required_columns} are required. Given dataframe has {columns_found}.".format(required_columns=str(required_columns), columns_found=str(list(any_df)))
        sys.exit("\n[QUIT]\n")
    # try to take other required columns, and fix if absent -> l
    try:
        l.append(any_df.loc[:,'shearsite'].to_frame())
    except:
        verbosePrint("[WARNING] checkEditDistance can't find 'shearsite' data. All the data will be processed as belonging to one unique shear site of nominal length '0'.")
        l.append(pd.DataFrame(data={'shearsite':['0']*len(any_df.loc[:,'randomBC'])}, index=any_df.loc[:,'randomBC'].index))
        # Note: here the index of 'fake shearsite data' must be forced to be the same of 'randomBC' due to the join in the following concat!!
        
    # Prepare DF: concat(l) --> DF
    DF = pd.concat(l, axis=1, join='inner')
    
    # check for correct usage - IS
    try:
        IS_set = set(any_df.loc[:,'genomic_coordinates'].tolist())
        if len(IS_set) != 1:
            verbosePrint("[WARNING] checkEditDistance found data from more than one IS found: {IS_list}. However, all the data will be processed as belonging to one unique IS!".format(IS_list=str(humanSorted(list(IS_set)))))
    except:
        verbosePrint("[WARNING] checkEditDistance cannot perform IS control on input data ('genomic_coordinates' not found). Thus, all the data will be processed as belonging to one unique IS!")
    # check for correct usage - Duplicate BCs (in shearsite groups, if given, or in general)
    a = DF.groupby('shearsite')['randomBC'].apply(lambda x: sorted(np.unique(x.tolist())))
    b = DF.groupby('shearsite')['randomBC'].apply(lambda x: sorted(x.tolist()))
    if not (a == b).all().all():
        verbosePrint("[WARNING] Duplicate randomBCs found (within shearsite groups, if given). Uniqueness will be forced through pandas.DataFrame.groupby('shearsite')['randomBC'].apply(lambda x: np.unique(x.tolist()))!")
    
    # Group randomBCs in lists, by shearsite, on DF
    DF = DF.groupby('shearsite')['randomBC'].apply(lambda x: np.unique(x.tolist()))
    # Here DF is returned as Series with Name=randomBC, Index=shearsite and values are lists of unique* randomBC (*nb: within the same list)
    # (np.unique is redundant in the use case for which this function was designed. Otherwise is forced (warnings if verbose!). See notes at the top
    
    # Collect edit distance matrix dataframes in l
    l = []
    l_append =l.append
    for shearsite, randomBC_list in DF.iteritems():
        # here kwargs 'groupLabel' are fundamental to l_append dataframes with mutually disjoint column labels
        # and then exploit (outer) 'join' at the end as a quick way to reshape data as a whole
        if all_combinations:
            inner_l = []
            inner_l_append = inner_l.append
            # compute all rows for current columns
            for shearsite2, randomBC_list2 in DF.iteritems():
                inner_l_append(editDistanceMatrix(randomBC_list, randomBC_list2, groupLabel1=buildGroupLabel(shearsite), groupLabel2=buildGroupLabel(shearsite2)))
            l_append(pd.concat(inner_l, axis=0, join='inner'))
        else:
            # compute only square sub-matrixes with rows equals to current columns
            l_append(editDistanceMatrix(randomBC_list, randomBC_list, groupLabel1=buildGroupLabel(shearsite), groupLabel2=buildGroupLabel(shearsite)))
    
    # Reshape data in a unique dataframe using multiple (outer) joins, exploiting unique columns, in any case:
    # editDistance_DF is a squared dataframe where rows and cols are equally labelled
    # (usually with randomBCs-seqs or with shearsite_randomBCs-seqs)
    editDistance_DF = pd.DataFrame().join(l, how='outer')
    return editDistance_DF

def editDistanceHeatmap(editDistance_DF, title='RANDOM-BARCODES EDIT-DISTANCE MATRIX', cmap="RdYlBu", annot=True, show_live=True, export=""):
    # Note: cmap = sns.diverging_palette(10, 133, l=60, n=12, center="dark", as_cmap=True) # red - dark - green
    # Note: with sns.palplot(sns.diverging_palette(10, 133, l=60, n=12, center="dark")) you can test and visualize palettes
    
    # get matrix dimension
    n = len(editDistance_DF.index)
    # Set up interactive mode (plot pop-upping)
    if show_live:
        plt.ion()
    else:
        plt.ioff()
    # Prepare plot environment
    min_size = 20
    max_size = 70
    x = int(ceil(min(max((n/2.0), min_size), max_size)))
    plt.figure(figsize=(x,x))
    plt.title(title)
    ax = plt.axes()
    # Control over annot kwarg
    annot_value = False
    annot_lim = 200
    if annot:
        if n>annot_lim:
            annot_value = False
            verbosePrint("[WARNING] editDistanceHeatmap cannot annotate matrixes with n>{annot_lim}.".format(annot_lim=str(annot_lim)))
        else:
            annot_value = True
    # Control for labels
    ticklabels = True
    alltick_lim = 400
    tick_lim = 800
    if n<tick_lim:
        if n>alltick_lim:
            ticklabels = int(ceil(n/100.0))
            verbosePrint("[WARNING] editDistanceHeatmap cannot write all labels of matrixes with n>{alltick_lim}. One each {ticklabels} will be reported instead.".format(alltick_lim=str(alltick_lim), ticklabels=str(ticklabels)))
    else:
        ticklabels = False
        verbosePrint("[WARNING] editDistanceHeatmap cannot write labels of matrixes with n>{tick_lim}. ".format(tick_lim=str(tick_lim)))
    # Set range for colors
    vmin=0
    vmax=12
    # Plot data
    sns.heatmap(editDistance_DF, ax=ax, cmap=cmap, annot=annot_value, xticklabels=ticklabels, yticklabels=ticklabels, vmin=vmin, vmax=vmax)  # ax = sns.heatmap(editDistance_DF, cmap=cmap)
    if show_live:
        plt.show()
    # Save plot
    if export:
        export = os.path.normpath(export)
        plt.savefig(export)
    # close
    if not show_live:
        plt.close()

def plotEditDistanceOccurrency(editDistance_DF, title="EDIT DISTANCE OCCURRENCIES", vmin=0, vmax=12, annot=True, percentile_colors=False, show_live=True, export=""):
    # Note: e.g. percentile_colors=[25, 75]
    # Note: vmin, vmax are the 'a-priori' min and max edit distances in editDistance_DF
    
    # Prepare data - NOTE: here data are DOUBLED because of the all-VS-all structure of editDistance_DF
    editDistance_array = np.array(editDistance_DF.values.tolist()).astype(float)
    np.fill_diagonal(editDistance_array, np.nan)  # set diagonal values to nan to remove uninformative 0's
    #np.arange(editDistance_array.shape[0])[:,None] > np.arange(editDistance_array.shape[1]) # mask: triangular upper
    editDistance_array[np.arange(editDistance_array.shape[0])[:,None] > np.arange(editDistance_array.shape[1])] = np.nan  # apply mask: triangular upper set to nan
    editDistance_array = editDistance_array[~np.isnan(editDistance_array)]  # remove nans
    editDistance_array = editDistance_array.astype(int)  # NOTE: here data are fixed: no doubled anymore, no uninformative zeros (diag) and non nan's
    # Set up interactive mode (plot pop-upping)
    if show_live:
        plt.ion()
    else:
        plt.ioff()
    # Prepare plot environment
    plt.figure(figsize=(13,13)) # Create matplotlib figure
    plt.title(title)
    ax = plt.axes() # Create matplotlib axes
    # Plot data
    counts, bins, patches = ax.hist(editDistance_array, facecolor='blue', edgecolor='black', bins=np.arange(start=vmin-0.5, stop=vmax+1.5, step=1))
    # Refine axis
    ax.set_ylabel("occurrencies (min={omin}, max={omax}, total={N})".format(omin=str(int(min(counts))), omax=str(int(max(counts))), N=str(int(sum(counts)))))
    ax.set_xlim([vmin-0.5,vmax+1.5])
    ax.set_xticks(range(vmin, vmax+1))
    # Optional: color bars according to percentiles
    if percentile_colors:
        # Change the colors of bars at the edges according to given percentiles!
        low_perc, high_perc = np.percentile(editDistance_array, percentile_colors)
        for patch, rightside, leftside in zip(patches, bins[1:], bins[:-1]):
            if rightside < low_perc:
                patch.set_facecolor('red')
            elif leftside > high_perc:
                patch.set_facecolor('green')
    # Optional: annotate x axis with counts and percentages
    if annot:
        # Label the raw counts and the percentages below the x-axis
        bin_centers = range(vmin, vmax+1)
        for count, x in zip(counts, bin_centers):
            # Label the raw counts
            count_label = str(int(count))
            ax.annotate(count_label, xy=(x, 0), xycoords=('data', 'axes fraction'), xytext=(0, -18), textcoords='offset points', va='top', ha='center')
            # Label the percentages
            percent = "{:.1f}".format(100 * float(count) / counts.sum()) + "%"
            ax.annotate(percent, xy=(x, 0), xycoords=('data', 'axes fraction'), xytext=(0, -32), textcoords='offset points', va='top', ha='center')
        # Refine axis
        ax.set_xlabel("edit distances\ncounts\npercentages")
        ax.xaxis.set_label_coords(1, -0.005)
        # Give more room at the bottom of the plot
        plt.subplots_adjust(bottom=0.2)
    else:
        # standard X axis label
        ax.set_xlabel("edit distances")
    # Show plot
    if show_live:
        plt.show()
    # Save plot
    if export:
        export = os.path.normpath(export)
        plt.savefig(export)
    # close
    if not show_live:
        plt.close()



def checkShearSitesOccurrency(any_df):
    # Note: here DISTINCT shearsites are counted through seq_count sum and/or by distinct randomBC count
    # data collector
    l = []
    # try to take required columns -> l
    required_columns = ['shearsite']
    try:
        l.append(any_df.loc[:,required_columns[0]].to_frame())
    except:
        print "\n[ERROR] checkShearSitesOccurrency wrong input! Column {required_column} is required. Given dataframe has {columns_found}.".format(required_column=str(required_columns[0]), columns_found=str(list(any_df)))
        sys.exit("\n[QUIT]\n")
    # try to take at least one 'quantification column'
    other_columns = ['seq_count', 'randomBC']
    columns_found = ['seq_count', 'randomBC']
    for c in other_columns:
        try:
            l.append(any_df.loc[:,c].to_frame())
        except:
            columns_found.remove(c)
    # Quit if input data are not sufficient
    if len(columns_found)<1:
        print "\n[ERROR] checkShearSitesOccurrency wrong input! At least one column among {other_columns} is required. Given dataframe has {columns_found}.".format(other_columns=str(other_columns), columns_found=str(list(any_df)))
        sys.exit("\n[QUIT]\n")
    # New DF with data to process
    DF = pd.concat(l, axis=1, join='inner')
    # results collector
    l = []
    if 'seq_count' in columns_found:
        l.append(pd.pivot_table(DF, columns='shearsite', values='seq_count', aggfunc=sum).to_frame())
    if 'randomBC' in columns_found:
        tmp_df = pd.pivot_table(DF, columns='shearsite', values='randomBC', aggfunc=lambda x: len(np.unique(x))).to_frame()
        tmp_df.rename(columns={'randomBC':'randomBC_count_distinct'}, inplace=True)
        l.append(tmp_df)
    # Build final data: dataframe with distinct shearsites as row index, seq_count and/or randomBC_count_distinct as columns, according to data availability
    ShearSitesOccurrency_DF = pd.concat(l, axis=1, join='inner')
    return ShearSitesOccurrency_DF

def plotShearSitesOccurrency(ShearSitesOccurrency_DF, title="SHEAR SITE OCCURRENCIES", normalize=100, show_live=True, export=""):
    # prepare data if required
    data = ShearSitesOccurrency_DF
    if normalize:
        data = ShearSitesOccurrency_DF.apply(lambda x: normalize*x/x.sum())
    # Set up interactive mode (plot pop-upping)
    if show_live:
        plt.ion()
    else:
        plt.ioff()
    # Prepare plot environment
    x = min(max(13, 13*len(ShearSitesOccurrency_DF)/50.0), 100)
    plt.figure(figsize=(x,13)) # Create matplotlib figure
    plt.title(title)
    ax = plt.axes() # Create matplotlib axes
    # Refine y axis
    # Refine axis
    tot_info = []
    if 'seq_count' in list(ShearSitesOccurrency_DF):
        tot_info.append("by SC={totBySC}".format(totBySC=str(ShearSitesOccurrency_DF['seq_count'].sum())))
    if 'randomBC_count_distinct' in list(ShearSitesOccurrency_DF):
        tot_info.append("by RBC countDistinct={totByRBCcd}".format(totByRBCcd=str(ShearSitesOccurrency_DF['randomBC_count_distinct'].sum())))
    tot_info = "[total: " + ", ".join(tot_info) + "]"
    if normalize:
        ax.set_ylabel("frequency (normalized to {normalize}) {tot_info}".format(normalize=str(normalize), tot_info=tot_info))
    else:
        ax.set_ylabel("occurrencies {tot_info}".format(tot_info=tot_info))
    # Plot data
    data.plot(ax=ax, kind='bar')
    # Show plot
    if show_live:
        plt.show()
    # Save plot
    if export:
        export = os.path.normpath(export)
        plt.savefig(export)
    # close
    if not show_live:
        plt.close()



def checkFragmentLengthDistribution(any_df):
    # Note: this function is designed to take in input data of many integration sites from a-sample/a-sonication-pool.
    #       HOWEVER, THIS FUNCTION ALWAYS WORKS, so be careful while choosing the input.
    # Note: this function compute the number of fragments in any_df (n of distinct couples [genomic_coordinates, shearsite]), ignoring any other
    #       column/hierarchy inside any_df
    
    # data collector
    l = []
    # try to take required columns -> l
    required_columns = ['shearsite', 'genomic_coordinates']
    try:
        [l.append(any_df.loc[:,c].to_frame()) for c in required_columns]
    except:
        print "\n[ERROR] checkFragmentLengthDistribution wrong input! Columns {required_columns} are required. Given dataframe has {columns_found}.".format(required_columns=str(required_columns), columns_found=str(list(any_df)))
        sys.exit("\n[QUIT]\n")
    # New DF with data to process
    DF = pd.concat(l, axis=1, join='inner')
    DF = DF.drop_duplicates()  #unique! -> locations_list, length_list = list(DF.genomic_coordinates), list(DF.shearsite)
    FragmentLengthDistribution_DF = pd.pivot_table(DF, columns='shearsite', values='genomic_coordinates', aggfunc=lambda x: len(x.unique())).to_frame()
    FragmentLengthDistribution_DF.index.name = required_columns[0]
    FragmentLengthDistribution_DF.rename(columns={required_columns[1]:'n_of_fragment'}, inplace=True)
    #FragmentLengthDistribution_DF.plot(kind='bar') # might be ok!
    #However, plotFragmentLengthDistribution is better!
    return FragmentLengthDistribution_DF
    
def plotFragmentLengthDistribution(FragmentLengthDistribution_DF, title= "FRAGMENT LENGTH DISTRIBUTION", binning='Freedman-Diaconis', normalize=True, rug=False, kde=True, sonicLengthEstimation=True, show_live=True, export=""):
    # Note about binning: None/False means 'bar plot' (no binning)
    #                     'Freedman-Diaconis' means 'autobinning'
    #                     Else, an int is expected
    # Note: some kwargs implies normalize=True overriding
    
    # override user input
    if kde is True:
        normalize = True
    if sonicLengthEstimation:
        normalize = True
    # reminpulate data
    exploded_data = []
    for ShS, n_of_fragment in FragmentLengthDistribution_DF.iterrows():
        exploded_data += [int(ShS)]*int(n_of_fragment)
    # Set up interactive mode (plot pop-upping)
    if show_live:
        plt.ion()
    else:
        plt.ioff()
    # Prepare plot environment
    plt.figure(figsize=(13,13)) # Create matplotlib figure
    plt.title(title)
    ax = plt.axes() # Create matplotlib axes
    # Add optional estimation of real distribution: sonicLengthEstimation
    # put here beacuse of legend issue
    if sonicLengthEstimation and (normalize or kde):
        import rpy2.robjects as robjects
        from rpy2.robjects.packages import importr
        sonicLength = importr("sonicLength")
        length_list = exploded_data  # shearsites exploded by count, i.e. the number of distinct genomic_coordinates
        locations_list = range(len(length_list))  # fake distinct genomic_coordinates
        results = sonicLength.estAbund(robjects.StrVector(locations_list), robjects.FloatVector(length_list))
        phi = results.rx2("phi")
        freq_phi = list(phi)
        length_phi = [float(l) for l in list(phi.names)]
        ax.plot(length_phi, freq_phi, 'r', lw=2, label="sonicLength\nestimated distribution")
    # Plot data
    bins=None
    if binning != 'Freedman-Diaconis':
        if (binning is False) or (binning is None):
            #bins = len(FragmentLengthDistribution_DF)
            lmin = min(FragmentLengthDistribution_DF.index.values.astype(int))
            lmax = max(FragmentLengthDistribution_DF.index.values.astype(int))
            bins = range(lmin, lmax+1)
        else:
            bins=binning
    ax = sns.distplot(exploded_data, bins=bins, norm_hist=normalize, rug=rug, kde=kde, kde_kws={'label': 'KDE - Gaussian kernel', "lw": 2}, hist_kws={'label': 'Observed Distribution', "color": "g"})
    # Refine axis
    if normalize:
        ax.set_ylabel("fragment length frequencies (total fragments={N}, occMin={omin}, occMax={omax})".format(omin=str(FragmentLengthDistribution_DF['n_of_fragment'].min()), omax=str(FragmentLengthDistribution_DF['n_of_fragment'].max()), N=str(FragmentLengthDistribution_DF['n_of_fragment'].sum())))
    else:
        ax.set_ylabel("fragment length occurrencies (min={omin}, max={omax}, total={N})".format(omin=str(FragmentLengthDistribution_DF['n_of_fragment'].min()), omax=str(FragmentLengthDistribution_DF['n_of_fragment'].max()), N=str(FragmentLengthDistribution_DF['n_of_fragment'].sum())))
    ax.set_xlabel("fragment lengths (ShearSites)")
    # Show plot
    if show_live:
        plt.show()
    # Save plot
    if export:
        export = os.path.normpath(export)
        plt.savefig(export)
    # close
    if not show_live:
        plt.close()
    
    

### IN DEVELOPMENT (?) #####################################################################
def checkShearSitesDistance(any_df):
    ShearSitesDistance_DF = None  # a matrix like the edit distance one, within a single is
    return ShearSitesDistance_DF
    
def plotShearSitesDistanceOccurrency(ShearSitesDistance_DF):
    # plot an hist like the one produced by plotEditDistanceOccurrency
    pass
############################################################################################



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
    # Up to now exhaustive_df is just like df with 'header_list' column added.
    # Keep in mind that this structure is ongoing and more columns will may appear. Let them be implicitly supported!!
    # Check advancements in matrix_RandomBC_processingModule.buildExhaustiveDataFrame

    
    #### Test Functions in this module ####
    show_live = True
    
    ## Barcodes nucleotide balancing
    how='distinct'
    overall_nucleotidesCount_DF = checkNucleotideBalancing(exhaustive_df, how=how, N_warn=True)  # or = checkNucleotideBalancing(df, ...)
    plotNucleotideBalancing(overall_nucleotidesCount_DF, title='[DEBUG] PILED-UP SEQUENCES'+" - "+how, show_live=show_live, export=os.path.join(os.getcwd(), "test_output", "debug_checkNucleotidesBalancing.pdf"))
#
#    ## Barcodes occurrencies
#    overall_distinctBC_DF = checkRandomBCoccurrency(exhaustive_df)
#    plotRandomBCoccurrency(overall_distinctBC_DF, title='[DEBUG] RANDOM-BARCODE OCCURRENCIES', show_live=show_live, export=os.path.join(os.getcwd(), "test_output", "debug_checkRandomBCoccurrency.pdf"))
#
#    ## Edit distance occurrencies
#    data = exhaustive_df.loc[:,'randomBC'].to_frame()
#    data = exhaustive_df.loc[:,['randomBC', 'shearsite']]
#    data = exhaustive_df
#    editDistance_DF = checkEditDistance(data, all_combinations=False)
#    plotEditDistanceOccurrency(editDistance_DF, title='[DEBUG] EDIT DISTANCE OCCURRENCIES', show_live=show_live, export=os.path.join(os.getcwd(), "test_output", "debug_checkEditDistance_plotEditDistanceOccurrency.pdf"))
#    
#    ## ShearSite occurrencies
#    data = exhaustive_df.loc[:50,['randomBC', 'shearsite']]
#    data = exhaustive_df.loc[:,['seq_count', 'shearsite']]
#    data = exhaustive_df.loc[:50,:]
#    ShearSitesOccurrency_DF = checkShearSitesOccurrency(data)
#    plotShearSitesOccurrency(ShearSitesOccurrency_DF, title='[DEBUG] SHEARSITE OCCURRENCIES', show_live=show_live, export=os.path.join(os.getcwd(), "test_output", "debug_checkShearSitesOccurrency.pdf"))
#
#    ## Fragment Length distribution plot
#    for barcode in exhaustive_df['barcode'].unique():
#        data = exhaustive_df[exhaustive_df['barcode'] == barcode]
#        FragmentLengthDistribution_DF = checkFragmentLengthDistribution(data)
#        plotFragmentLengthDistribution(FragmentLengthDistribution_DF, title='[DEBUG] FRAGMENT LENGTH DISTRIBUTION', show_live=show_live, export=os.path.join(os.getcwd(), "test_output", "debug_{barcode}_FragmentLengthDistribution.pdf".format(barcode=barcode)))
#
    ## Edit distance matrixes

#    all_combinations=False
#    data = exhaustive_df.loc[:,['randomBC', 'shearsite']]
#    editDistance_DF = checkEditDistance(data.loc[:0,:], all_combinations=all_combinations)
#    editDistanceHeatmap(editDistance_DF, title='RANDOM-BARCODES EDIT-DISTANCE MATRIX', cmap="RdYlBu", annot=True, show_live=False, export=os.path.join(os.getcwd(), "test_output", "debug_checkEditDistance_diagonal_othercolors_0.pdf"))
#    editDistance_DF = checkEditDistance(data.loc[:2,:], all_combinations=all_combinations)
#    editDistanceHeatmap(editDistance_DF, title='RANDOM-BARCODES EDIT-DISTANCE MATRIX', cmap="RdYlBu", annot=True, show_live=False, export=os.path.join(os.getcwd(), "test_output", "debug_checkEditDistance_diagonal_othercolors_2.pdf"))
#    editDistance_DF = checkEditDistance(data.loc[:5,:], all_combinations=all_combinations)
#    editDistanceHeatmap(editDistance_DF, title='RANDOM-BARCODES EDIT-DISTANCE MATRIX', cmap="RdYlBu", annot=True, show_live=False, export=os.path.join(os.getcwd(), "test_output", "debug_checkEditDistance_diagonal_othercolors_5.pdf"))
#    editDistance_DF = checkEditDistance(data.loc[:10,:], all_combinations=all_combinations)
#    editDistanceHeatmap(editDistance_DF, title='RANDOM-BARCODES EDIT-DISTANCE MATRIX', cmap="RdYlBu", annot=True, show_live=False, export=os.path.join(os.getcwd(), "test_output", "debug_checkEditDistance_diagonal_othercolors_10.pdf"))
#    editDistance_DF = checkEditDistance(data.loc[:25,:], all_combinations=all_combinations)
#    editDistanceHeatmap(editDistance_DF, title='RANDOM-BARCODES EDIT-DISTANCE MATRIX', cmap="RdYlBu", annot=True, show_live=False, export=os.path.join(os.getcwd(), "test_output", "debug_checkEditDistance_diagonal_othercolors_25.pdf"))
#    editDistance_DF = checkEditDistance(data.loc[:50,:], all_combinations=all_combinations)
#    editDistanceHeatmap(editDistance_DF, title='RANDOM-BARCODES EDIT-DISTANCE MATRIX', cmap="RdYlBu", annot=True, show_live=False, export=os.path.join(os.getcwd(), "test_output", "debug_checkEditDistance_diagonal_othercolors_50.pdf"))
#    editDistance_DF = checkEditDistance(data.loc[:75,:], all_combinations=all_combinations)
#    editDistanceHeatmap(editDistance_DF, title='RANDOM-BARCODES EDIT-DISTANCE MATRIX', cmap="RdYlBu", annot=True, show_live=False, export=os.path.join(os.getcwd(), "test_output", "debug_checkEditDistance_diagonal_othercolors_75.pdf"))
#    editDistance_DF = checkEditDistance(data.loc[:100,:], all_combinations=all_combinations)
#    editDistanceHeatmap(editDistance_DF, title='RANDOM-BARCODES EDIT-DISTANCE MATRIX', cmap="RdYlBu", annot=True, show_live=False, export=os.path.join(os.getcwd(), "test_output", "debug_checkEditDistance_diagonal_othercolors_100.pdf"))
#    editDistance_DF = checkEditDistance(data.loc[:150,:], all_combinations=all_combinations)
#    editDistanceHeatmap(editDistance_DF, title='RANDOM-BARCODES EDIT-DISTANCE MATRIX', cmap="RdYlBu", annot=True, show_live=False, export=os.path.join(os.getcwd(), "test_output", "debug_checkEditDistance_diagonal_othercolors_150.pdf"))
#    editDistance_DF = checkEditDistance(data.loc[:200,:], all_combinations=all_combinations)
#    editDistanceHeatmap(editDistance_DF, title='RANDOM-BARCODES EDIT-DISTANCE MATRIX', cmap="RdYlBu", annot=True, show_live=False, export=os.path.join(os.getcwd(), "test_output", "debug_checkEditDistance_diagonal_othercolors_200.pdf"))
#    editDistance_DF = checkEditDistance(data.loc[:300,:], all_combinations=all_combinations)
#    editDistanceHeatmap(editDistance_DF, title='RANDOM-BARCODES EDIT-DISTANCE MATRIX', cmap="RdYlBu", annot=True, show_live=False, export=os.path.join(os.getcwd(), "test_output", "debug_checkEditDistance_diagonal_othercolors_300.pdf"))
#    editDistance_DF = checkEditDistance(data.loc[:400,:], all_combinations=all_combinations)
#    editDistanceHeatmap(editDistance_DF, title='RANDOM-BARCODES EDIT-DISTANCE MATRIX', cmap="RdYlBu", annot=True, show_live=False, export=os.path.join(os.getcwd(), "test_output", "debug_checkEditDistance_diagonal_othercolors_400.pdf"))
#    editDistance_DF = checkEditDistance(data.loc[:450,:], all_combinations=all_combinations)
#    editDistanceHeatmap(editDistance_DF, title='RANDOM-BARCODES EDIT-DISTANCE MATRIX', cmap="RdYlBu", annot=True, show_live=False, export=os.path.join(os.getcwd(), "test_output", "debug_checkEditDistance_diagonal_othercolors_450.pdf"))
#    editDistance_DF = checkEditDistance(data.loc[:500,:], all_combinations=all_combinations)
#    editDistanceHeatmap(editDistance_DF, title='RANDOM-BARCODES EDIT-DISTANCE MATRIX', cmap="RdYlBu", annot=True, show_live=False, export=os.path.join(os.getcwd(), "test_output", "debug_checkEditDistance_diagonal_othercolors_500.pdf"))
#    editDistance_DF = checkEditDistance(data.loc[:550,:], all_combinations=all_combinations)
#    editDistanceHeatmap(editDistance_DF, title='RANDOM-BARCODES EDIT-DISTANCE MATRIX', cmap="RdYlBu", annot=True, show_live=False, export=os.path.join(os.getcwd(), "test_output", "debug_checkEditDistance_diagonal_othercolors_550.pdf"))
#    editDistance_DF = checkEditDistance(data.loc[:600,:], all_combinations=all_combinations)
#    editDistanceHeatmap(editDistance_DF, title='RANDOM-BARCODES EDIT-DISTANCE MATRIX', cmap="RdYlBu", annot=True, show_live=False, export=os.path.join(os.getcwd(), "test_output", "debug_checkEditDistance_diagonal_othercolors_600.pdf"))
#    editDistance_DF = checkEditDistance(data, all_combinations=all_combinations)
#    editDistanceHeatmap(editDistance_DF, title='RANDOM-BARCODES EDIT-DISTANCE MATRIX', cmap="RdYlBu", annot=True, show_live=False, export=os.path.join(os.getcwd(), "test_output", "debug_checkEditDistance_diagonal_othercolors_all.pdf"))
#
#    all_combinations=True
#    data = exhaustive_df.loc[:,['randomBC', 'shearsite']]
#    editDistance_DF = checkEditDistance(data.loc[:0,:], all_combinations=all_combinations)
#    editDistanceHeatmap(editDistance_DF, title='RANDOM-BARCODES EDIT-DISTANCE MATRIX', cmap="RdYlBu", annot=True, show_live=False, export=os.path.join(os.getcwd(), "test_output", "debug_checkEditDistance_othercolors_0.pdf"))
#    editDistance_DF = checkEditDistance(data.loc[:2,:], all_combinations=all_combinations)
#    editDistanceHeatmap(editDistance_DF, title='RANDOM-BARCODES EDIT-DISTANCE MATRIX', cmap="RdYlBu", annot=True, show_live=False, export=os.path.join(os.getcwd(), "test_output", "debug_checkEditDistance_othercolors_2.pdf"))
#    editDistance_DF = checkEditDistance(data.loc[:5,:], all_combinations=all_combinations)
#    editDistanceHeatmap(editDistance_DF, title='RANDOM-BARCODES EDIT-DISTANCE MATRIX', cmap="RdYlBu", annot=True, show_live=False, export=os.path.join(os.getcwd(), "test_output", "debug_checkEditDistance_othercolors_5.pdf"))
#    editDistance_DF = checkEditDistance(data.loc[:10,:], all_combinations=all_combinations)
#    editDistanceHeatmap(editDistance_DF, title='RANDOM-BARCODES EDIT-DISTANCE MATRIX', cmap="RdYlBu", annot=True, show_live=False, export=os.path.join(os.getcwd(), "test_output", "debug_checkEditDistance_othercolors_10.pdf"))
#    editDistance_DF = checkEditDistance(data.loc[:25,:], all_combinations=all_combinations)
#    editDistanceHeatmap(editDistance_DF, title='RANDOM-BARCODES EDIT-DISTANCE MATRIX', cmap="RdYlBu", annot=True, show_live=False, export=os.path.join(os.getcwd(), "test_output", "debug_checkEditDistance_othercolors_25.pdf"))
#    editDistance_DF = checkEditDistance(data.loc[:50,:], all_combinations=all_combinations)
#    editDistanceHeatmap(editDistance_DF, title='RANDOM-BARCODES EDIT-DISTANCE MATRIX', cmap="RdYlBu", annot=True, show_live=False, export=os.path.join(os.getcwd(), "test_output", "debug_checkEditDistance_othercolors_50.pdf"))
#    editDistance_DF = checkEditDistance(data.loc[:75,:], all_combinations=all_combinations)
#    editDistanceHeatmap(editDistance_DF, title='RANDOM-BARCODES EDIT-DISTANCE MATRIX', cmap="RdYlBu", annot=True, show_live=False, export=os.path.join(os.getcwd(), "test_output", "debug_checkEditDistance_othercolors_75.pdf"))
#    editDistance_DF = checkEditDistance(data.loc[:100,:], all_combinations=all_combinations)
#    editDistanceHeatmap(editDistance_DF, title='RANDOM-BARCODES EDIT-DISTANCE MATRIX', cmap="RdYlBu", annot=True, show_live=False, export=os.path.join(os.getcwd(), "test_output", "debug_checkEditDistance_othercolors_100.pdf"))
#    editDistance_DF = checkEditDistance(data.loc[:150,:], all_combinations=all_combinations)
#    editDistanceHeatmap(editDistance_DF, title='RANDOM-BARCODES EDIT-DISTANCE MATRIX', cmap="RdYlBu", annot=True, show_live=False, export=os.path.join(os.getcwd(), "test_output", "debug_checkEditDistance_othercolors_150.pdf"))
#    editDistance_DF = checkEditDistance(data.loc[:200,:], all_combinations=all_combinations)
#    editDistanceHeatmap(editDistance_DF, title='RANDOM-BARCODES EDIT-DISTANCE MATRIX', cmap="RdYlBu", annot=True, show_live=False, export=os.path.join(os.getcwd(), "test_output", "debug_checkEditDistance_othercolors_200.pdf"))
#    editDistance_DF = checkEditDistance(data.loc[:300,:], all_combinations=all_combinations)
#    editDistanceHeatmap(editDistance_DF, title='RANDOM-BARCODES EDIT-DISTANCE MATRIX', cmap="RdYlBu", annot=True, show_live=False, export=os.path.join(os.getcwd(), "test_output", "debug_checkEditDistance_othercolors_300.pdf"))
#    editDistance_DF = checkEditDistance(data.loc[:400,:], all_combinations=all_combinations)
#    editDistanceHeatmap(editDistance_DF, title='RANDOM-BARCODES EDIT-DISTANCE MATRIX', cmap="RdYlBu", annot=True, show_live=False, export=os.path.join(os.getcwd(), "test_output", "debug_checkEditDistance_othercolors_400.pdf"))
#    editDistance_DF = checkEditDistance(data.loc[:450,:], all_combinations=all_combinations)
#    editDistanceHeatmap(editDistance_DF, title='RANDOM-BARCODES EDIT-DISTANCE MATRIX', cmap="RdYlBu", annot=True, show_live=False, export=os.path.join(os.getcwd(), "test_output", "debug_checkEditDistance_othercolors_450.pdf"))
#    editDistance_DF = checkEditDistance(data.loc[:500,:], all_combinations=all_combinations)
#    editDistanceHeatmap(editDistance_DF, title='RANDOM-BARCODES EDIT-DISTANCE MATRIX', cmap="RdYlBu", annot=True, show_live=False, export=os.path.join(os.getcwd(), "test_output", "debug_checkEditDistance_othercolors_500.pdf"))
#    editDistance_DF = checkEditDistance(data.loc[:550,:], all_combinations=all_combinations)
#    editDistanceHeatmap(editDistance_DF, title='RANDOM-BARCODES EDIT-DISTANCE MATRIX', cmap="RdYlBu", annot=True, show_live=False, export=os.path.join(os.getcwd(), "test_output", "debug_checkEditDistance_othercolors_550.pdf"))
#    editDistance_DF = checkEditDistance(data.loc[:600,:], all_combinations=all_combinations)
#    editDistanceHeatmap(editDistance_DF, title='RANDOM-BARCODES EDIT-DISTANCE MATRIX', cmap="RdYlBu", annot=True, show_live=False, export=os.path.join(os.getcwd(), "test_output", "debug_checkEditDistance_othercolors_600.pdf"))
#    editDistance_DF = checkEditDistance(data, all_combinations=all_combinations)
#    editDistanceHeatmap(editDistance_DF, title='RANDOM-BARCODES EDIT-DISTANCE MATRIX', cmap="RdYlBu", annot=True, show_live=False, export=os.path.join(os.getcwd(), "test_output", "debug_checkEditDistance_othercolors_all.pdf"))
#
#    # Any df is accepted so the contents of each plot is the whole input DF.
#    #title and export kwargs allows you to put the proper label to your contents!