# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19 10:11:30 2016

@author: stefano
"""

###################################################################################
### HERE A 'DATASET' IS SUPPOSED TO BE UNIVOCALLY DETERMINED BY (DISEASE, PATIENT)
### AND TO CONTAIN AT LEAST ONE POOL
###################################################################################

###################################################################################
### HERE A 'DATASET TUPLE' IS SUPPOSED TO BE:
### (ground_dir, DISEASE, PATIENT, False) - all pools available
### or
### (ground_dir, DISEASE, PATIENT, specific_POOLS list)
### Thus data and paths result fully specified
###################################################################################



#++++++++++++++++++++++++++++ Requested Package(s) Import +++++++++++++++++++++++++++++#
import os, sys
import pprint
import matrix_configure_module


#++++++++++++++++++++++ Global Vars +++++++++++++++++++++++#
verbose = matrix_configure_module.verbose


#++++++++++++++++++++++ Global Funcs ++++++++++++++++++++++#
verbosePrint = matrix_configure_module.verbosePrint
humanSorted = matrix_configure_module.humanSorted


#+++++++++++++++++++++++++++++++++++++++ FUNCTIONS +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#


def getBASEDIR(ground_dir, DISEASE, PATIENT):
    '''
    PURPOSE
    Build and return the standard BASEDIR path, regularly built as:
    $ground_dir/$DISEASE/$PATIENT/
    checked for existence and readability
    '''
    # Check ground_dir
    ground_dir = str(ground_dir)
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
    # Build BASEDIR
    DISEASE = str(DISEASE)
    PATIENT = str(PATIENT)
    try:
        BASEDIR = os.path.normpath(os.path.join(ground_dir, DISEASE, PATIENT))
    except Exception, err_message:
        print "\n[ERROR] Unable to build BASEDIR as a valid path!"
        print "ground_dir='{ground_dir}'".format(ground_dir=str(ground_dir))
        print "DISEASE='{DISEASE}'".format(DISEASE=str(DISEASE))
        print "PATIENT='{PATIENT}'".format(PATIENT=str(PATIENT))
        print "os.path.normpath(os.path.join(ground_dir, DISEASE, PATIENT)) returned: ", err_message
        sys.exit("\n[QUIT]\n")
    # Check BASEDIR
    if not os.path.exists(BASEDIR):
        print "\n[ERROR] BASEDIR='{BASEDIR}' does not exist!".format(BASEDIR=str(BASEDIR))
        sys.exit("\n[QUIT]\n")
    if not os.access(BASEDIR, os.R_OK):
        print "\n[ERROR] You have not read permission in BASEDIR='{BASEDIR}'".format(BASEDIR=str(BASEDIR))
        sys.exit("\n[QUIT]\n")
    return BASEDIR
    

def getContentdir(BASEDIR, content):
    '''
    PURPOSE
    Build and return the standard 'content' path, regularly built as:
    $BASEDIR/$content
    checked for existence and readability
    
    NOTE: 'content' is something like: bam, bed, iss, quality, quantification
    '''
    # Build contentdir
    content = str(content)
    try:
        contentdir = os.path.normpath(os.path.join(BASEDIR, content))
    except Exception, err_message:
        print "\n[ERROR] Unable to build contentdir as a valid path!"
        print "BASEDIR='{BASEDIR}'".format(BASEDIR=str(BASEDIR))
        print "content='{content}'".format(content=str(content))
        print "os.path.normpath(os.path.join(BASEDIR, content)) returned: ", err_message
        sys.exit("\n[QUIT]\n")
    # Check contentdir
    if not os.path.exists(contentdir):
        print "\n[ERROR] contentdir='{contentdir}' does not exist!".format(contentdir=str(contentdir))
        sys.exit("\n[QUIT]\n")
    if not os.access(contentdir, os.R_OK):
        print "\n[ERROR] You have not read permission in contentdir='{contentdir}'".format(contentdir=str(contentdir))
        sys.exit("\n[QUIT]\n")
    return contentdir


def getPOOLlist(contentdir):
    '''
    PURPOSE
    Return pool list as list of folder names in contentdir
    '''
    contentdir_ls = os.listdir(contentdir)
    POOL_list = [POOL for POOL in contentdir_ls if os.path.isdir(os.path.join(contentdir, POOL))]
    POOL_list = humanSorted(POOL_list)
    return POOL_list


def listDir(folder, name_filter=None):
    """
    PURPOSE
    Returns the content (abs path list) of input folder, filtered
    by name_filter substring, if not None.
    IN:
    folder - ... complete path (e.g. starting with '/') are checked and taken as they are.
                 words or partial path ('foldername' or 'foldername/subfolder') will be
                 taken as cwd/...
    name_filter - substring for exact pattern matching / None for no filtering
    OUT:
    dir_content or filtered_dir_content
    """
    # Check folder
    try:
        folder = os.path.abspath(os.path.normpath(folder))
    except Exception, err_message:
        print "\n[ERROR] folder='{folder}' is not formatted as valid path!".format(folder=str(folder))
        print "os.path.normpath returned: ", err_message
        sys.exit("\n[QUIT]\n")
    if os.path.isfile(folder):
        print "\n[ERROR] folder must be a folder path, not a file! Your input: folder='{folder}'".format(folder=str(folder))
        sys.exit("\n[QUIT]\n")
    if not os.path.exists(folder):
        print "\n[ERROR] folder='{folder}' does not exist!".format(folder=str(folder))
        sys.exit("\n[QUIT]\n")
    if not os.access(folder, os.R_OK):
        print "\n[ERROR] You have not read permission in folder='{folder}'".format(folder=str(folder))
        sys.exit("\n[QUIT]\n")
    # Get folder contents
    dir_content = []
    for filename in os.listdir(folder):
        dir_content.append(os.path.normpath(os.path.join(folder, filename)))
    if name_filter is None:
        return humanSorted(dir_content)
    else:
        filtered_dir_content = []
        for line in dir_content:
            if name_filter in line:
                filtered_dir_content.append(line)
        return humanSorted(filtered_dir_content)


def getDatasetPathDict (ground_dir, DISEASE, PATIENT, specific_POOLS=False):
    '''
    PURPOSE
    Return dataset_path_dict = {'DISEASE': {*}}
                         {*} = {'PATIENT': {**}}
                        {**} = {'POOL': {***}}
                       {***} = {'source_data_dir': absolute path to refacored-tsv-data FOLDER,
                                'randomBC_data_dir': absolute path to randomBC-fasta FOLDER,
                                'source_data_path': absolute path of refacored-tsv-data FILE, (A)
                                'randomBC_data_path': absolute path to randomBC-fasta FILE,   (A)
                                }
    INPUT NOTE: 
    - specific_POOLS must be False (take all available pools) or a LIST of existing pools
    OUTPUT NOTE:
    (A) returned only if found!
    '''
    
    def BuildPathAndCheck(abs_path, POOL):
        try:
            POOLpath = os.path.normpath(os.path.join(abs_path, str(POOL)))
        except Exception, err_message:
            print "\n[ERROR] Unable to build a valid POOLpath for POOL='{POOL}'!".format(POOL=str(POOL))
            print "abs_path='{abs_path}'".format(abs_path=str(abs_path))
            print "os.path.normpath(os.path.join(abs_path, str(POOL))) returned: ", err_message
            sys.exit("\n[QUIT]\n")
        if not os.path.exists(POOLpath):
            print "\n[ERROR] POOLpath='{POOLpath}' does not exist!".format(POOLpath=str(POOLpath))
            sys.exit("\n[QUIT]\n")
        if not os.access(POOLpath, os.R_OK):
            print "\n[ERROR] You have not read permission in POOLpath='{POOLpath}'".format(POOLpath=str(POOLpath))
            sys.exit("\n[QUIT]\n")
        return POOLpath
        
    # Build BASEDIR and check
    BASEDIR = getBASEDIR(ground_dir, DISEASE, PATIENT)
    # Build contentdirs and check
    refactored_data_dir = getContentdir(BASEDIR, 'iss')
    randomBC_data_dir = getContentdir(BASEDIR, 'quality')

    # Get POOL_list (explicitly, as specific_POOLS kwarg, otherwise all pools found in 'refactored_data_dir' will be taken)
    POOL_list = None
    if specific_POOLS is False:
        # All pools
        POOL_list = getPOOLlist(refactored_data_dir)
    else:
        # Some pools
        if type(specific_POOLS) is not list:
            print "\n[ERROR] specific_POOLS must be a LIST or set to False. Your input: '{specific_POOLS}', type={t}".format(specific_POOLS=str(specific_POOLS), t=str(type(specific_POOLS)))
            sys.exit("\n[QUIT]\n")
        if (set(getPOOLlist(refactored_data_dir)).intersection(set(specific_POOLS)) != set(specific_POOLS)):
            print "\n[ERROR] specific_POOLS must be a list of EXISTING POOLS!"
            print "          * POOLS found in '{refactored_data_dir}': ".format(refactored_data_dir=str(refactored_data_dir)), humanSorted(getPOOLlist(refactored_data_dir))
            print "          * Yuor input POOLS: ", humanSorted(specific_POOLS)
            sys.exit("\n[QUIT]\n")
        POOL_list = specific_POOLS
    if not POOL_list:
        print "\n[ERROR] Can't find any POOL in refactored_data_dir='{refactored_data_dir}'".format(refactored_data_dir=str(refactored_data_dir))
        sys.exit("\n[QUIT]\n")
    
    # Fill dataset_path_dict with 'source_data_dir' and 'randomBC_data_dir'
    dataset_path_dict = {DISEASE:{PATIENT:{}}}
    for POOL in POOL_list:
        dataset_path_dict[DISEASE][PATIENT][POOL] = {}
        dataset_path_dict[DISEASE][PATIENT][POOL]['source_data_dir'] = BuildPathAndCheck(refactored_data_dir, POOL)
        dataset_path_dict[DISEASE][PATIENT][POOL]['randomBC_data_dir'] = BuildPathAndCheck(randomBC_data_dir, POOL)
    
        # Looking into source_data_dir for source_data_path and fill dataset_path_dict
        l = listDir(dataset_path_dict[DISEASE][PATIENT][POOL]['source_data_dir'], name_filter='_refactored.tsv.gz')
        if len(l) != 1:
            print "\n[ERROR] Unexpected content of dir", dataset_path_dict[DISEASE][PATIENT][POOL]['source_data_dir']
            sys.exit("\n[QUIT]\n")
        else:
            dataset_path_dict[DISEASE][PATIENT][POOL]['source_data_path'] = os.path.normpath(l[0])
        
        # Looking into randomBC_data_dir for randomBC_data_path  and fill dataset_path_dict
        l = listDir(dataset_path_dict[DISEASE][PATIENT][POOL]['randomBC_data_dir'], name_filter='.fa.gz')
        if len(l) != 1:
            if len(l) == 0:
                # randomBC not available
                verbosePrint("[Warning] Random Barcodes not available for '{DISEASE}-{PATIENT}-{POOL}'!".format(DISEASE=str(DISEASE), PATIENT=str(PATIENT), POOL=str(POOL)))
                del dataset_path_dict[DISEASE][PATIENT][POOL]['randomBC_data_dir']
            else:
                print "\n[ERROR] Unexpected content of dir", dataset_path_dict[DISEASE][PATIENT][POOL]['randomBC_data_dir']
                sys.exit("\n[QUIT]\n")
        else:
            dataset_path_dict[DISEASE][PATIENT][POOL]['randomBC_data_path'] = os.path.normpath(l[0])

    return dataset_path_dict


def getLaunchPathDict(dataset_tuple_list):
    '''
    PURPOSE
    Take a list of dataset tuple, loop and call getDatasetPathDict,
    join returned dataset_path_dicts and finally return a comprehensive
    dict for the whole launch, called launch_path_dict.
    
    NOTE:
    dataset_tuple_list is like:
    [(a_ground_dir, somewhat_DISEASE, a_PATIENT, False), (a_ground_dir, somewhat_DISEASE, another_PATIENT, a_list_of_pools_of_interest), ... so on ...]
    '''
    verbosePrint("\n>>> Searching for data sources...")

    launch_path_dict = {}
    for ground_dir, DISEASE, PATIENT, specific_POOLS in dataset_tuple_list:
        dataset_path_dict = getDatasetPathDict(ground_dir, DISEASE, PATIENT, specific_POOLS)
        if DISEASE not in launch_path_dict:
            launch_path_dict[DISEASE] = dataset_path_dict[DISEASE]
        else:
            if PATIENT not in launch_path_dict[DISEASE]:
                launch_path_dict[DISEASE][PATIENT] = dataset_path_dict[DISEASE][PATIENT]
            else:
                print "\n[ERROR] Duplicate DISEASE-PATIENT found!", (DISEASE, PATIENT)
                print "Your input dataset_tuple_list:", dataset_tuple_list
                sys.exit("\n[QUIT]\n")
    
    verbosePrint("\n> Data found:")
    if verbose:
        pprint.pprint(launch_path_dict)
    verbosePrint(">>> Done!")

    return launch_path_dict


#++++++++++++++++++++++++++++++++++++++ MAIN and TEST +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#


if __name__ == "__main__":
    
    ### Test vars ###
    # 2 dataset:
    #   - Ferrari, thalassemia, all pools
    #   - HIV, circulatingDNA, pool 'cLR1' only
    dataset_tuple_list = [ ('/opt/NGS/results', 'Ferrari', 'PE_Thalassemia_inVitro', False), ('/opt/NGS/results', 'HIV_patients', 'circDNA', ['cLR1']) ]
    
    ### Test code ###
    launch_path_dict = getLaunchPathDict(dataset_tuple_list)

