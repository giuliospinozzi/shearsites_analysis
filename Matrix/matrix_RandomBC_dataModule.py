# -*- coding: utf-8 -*-
"""
Created on Fri Jan  8 10:20:07 2016

@author: stefano
"""

#++++++++++++++ Requested Package(s) Import +++++++++++++++#
import os, sys, csv
import matrix_RandomBC_globModule


#++++++++++++++++++++++ Global Vars +++++++++++++++++++++++#
verbose = matrix_RandomBC_globModule.verbose


#++++++++++++++++++++++ Global Funcs ++++++++++++++++++++++#
verbosePrint = matrix_RandomBC_globModule.verbosePrint
humanSorted = matrix_RandomBC_globModule.humanSorted


#+++++++++++++++++++++++++++++++++++++++ FUNCTIONS +++++++++++++++++++++++++++++++++++++++#

def buildInputPath(ground_dir, DISEASE, PATIENT, POOL):
    '''
    Purpose:
    Build and return the standard input path (thought as containing all the input file<->sample<->barcode of interest),
    starting from ground_dir
    IN:
    ground_dir (abs valid path), DISEASE, PATIENT, POOL (strings)
    OUT: INDIR = $ground_dir/$DISEASE/$PATIENT/'quantification'/$POOL/'RandomBC' (checked)
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
    INDIR = os.path.normpath(os.path.join(BASEDIR, "quantification", POOL, "RandomBC"))
    # Check INDIR
    if not os.path.exists(INDIR):
        print "\n[ERROR] INDIR='{INDIR}' just build does not exist!".format(INDIR=str(INDIR))
        sys.exit("\n[QUIT]\n")
    if not os.access(INDIR, os.R_OK):
        print "\n[ERROR] You have not read permission in INDIR='{INDIR}'".format(INDIR=str(INDIR))
        sys.exit("\n[QUIT]\n")
    return INDIR


def listDir(folder, name_filter=None):
    """
    Purpose:
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


def loadFile (path, delimiter):
    """
    Purpose: LOAD A delimiter-DELIMITED FILE AS NESTED LIST
    IN:
    path, delimiter
    OUT:
    nested_list - nested list of file contents (row > columns)
    """
    nested_list = []
    nested_list_append = nested_list.append   # improves performaces
    with open(path) as csvfile:
        reader = csv.reader(csvfile, delimiter=delimiter)
        for row in reader:
            nested_list_append(row)
    return nested_list
    
    
def arrangeData(data_file_nested_list):
    """
    Purpose: take in input data as nested-list from loadFile function and return two dict as shown below;
             Data in nested-list are supposed to come from a single/specific file<->sample<->barcode!
    IN:
    data_file_nested_list (as returned by loadFile, for a single/specific file<->sample<->barcode)
    OUT:
    alldata_dict = {'IS_key': line_dict, ...}
                              line_dict = {'header': line_data, ...}
                                                     line_data = {'r1_chr': line[1].replace("chr", ""),
                                                                  'r1_start': int(line[2]),
                                                                  'r1_end': int(line[3]),
                                                                  'r1_quality': int(line[4]),
                                                                  'r1_strand': line[5],
                                                                  'r2_chr': line[6].replace("chr", ""),
                                                                  'r2_start': int(line[7]),
                                                                  'r2_end': int(line[8]),
                                                                  'r2_quality': int(line[9]),
                                                                  'r2_strand': line[10],
                                                                  'length': int(line[11]),
                                                                  'randomBC': line[12],
                                                                  }
                                                                  
    IS_dict = {'IS_key': shearsite_dict, ...}
                         shearsite_dict = {'shearsite_key': randomBC_dict, ...}
                                                            randomBC_dict = {'randomBC_key': barcode_count, ...}
    NOTE: - alldata_dict structure is designed to provide the input for matrix_RandomBC_processingModule.buildExhaustiveDataFrame({'barcode': alldata_dict, ...})
            in order to create exhaustive_df
          - IS_dict structure is designed to provide the input for matrix_RandomBC_processingModule.buildDataFrame({'barcode': IS_dict, ...})
            in order to create df 
    """
    def generateISkey(stuff_list, join_char='_'):
        stuff_list = [str(x) for x in stuff_list]
        return join_char.join(stuff_list)
    def generateShearsiteKey(item):
        return item
    def generateRandomBCkey(item):
        return item
    
    IS_dict = {}
    alldata_dict = {}
    for line in data_file_nested_list:
        
        # Load all data
        line_dict = {}
        header = line[0]
        line_dict[header] = {'r1_chr': line[1].replace("chr", ""),
                             'r1_start': int(line[2]),
                             'r1_end': int(line[3]),
                             'r1_quality': int(line[4]),
                             'r1_strand': line[5],
                             'r2_chr': line[6].replace("chr", ""),
                             'r2_start': int(line[7]),
                             'r2_end': int(line[8]),
                             'r2_quality': int(line[9]),
                             'r2_strand': line[10],
                             'length': int(line[11]),
                             'randomBC': line[12],
                            }
        IS_key = generateISkey([line_dict[header]['r1_chr'], line_dict[header]['r1_start'], line_dict[header]['r1_strand']])
        if alldata_dict.get(IS_key) is None:
            alldata_dict[IS_key] = line_dict
        else:
            alldata_dict[IS_key].update(line_dict)
        
        # IS layer
        #IS_key = generateISkey()  # already done above
        shearsite_dict = None
        
        if IS_dict.get(IS_key) is None:
            IS_dict[IS_key] = {}
        shearsite_dict = IS_dict[IS_key]

        # ShearSite layer
        shearsite_key = generateShearsiteKey(line_dict[header]['length'])
        randomBC_dict = None
        
        if shearsite_dict.get(shearsite_key) is None:
            shearsite_dict[shearsite_key] = {}
        randomBC_dict = shearsite_dict[shearsite_key]
        
        # RandomBC layer
        randomBC_key = generateRandomBCkey(line_dict[header]['randomBC'])
        if randomBC_dict.get(randomBC_key) is None:
            randomBC_dict[randomBC_key] = 1
        else:
            randomBC_dict[randomBC_key] += 1
    
    return alldata_dict, IS_dict
        

def loadDataFiles(ground_dir, DISEASE, PATIENT, POOL, data_files_name_filter, data_files_delimiter):
    """
    Purpose: Implement loadFile and arrangeData in a loop over all data files (all barcodes) of the POOL.
             Results are returned in two dict: POOL_alldata_dict, POOL_IS_dict; both have barcodes as keys.
    """
    verbosePrint("\n>>> Loading data ...")
    INDIR = buildInputPath(ground_dir, DISEASE, PATIENT, POOL)
    verbosePrint("> path: {path}".format(path=str(INDIR)))
    filtered_dir_content = listDir(INDIR, name_filter=data_files_name_filter)
    verbosePrint("> exploited substring for data detection: '{data_files_name_filter}'".format(data_files_name_filter=str(data_files_name_filter)))
    verbosePrint("> n data files detected: {n_files}".format(n_files=str(len(filtered_dir_content))))
    verbosePrint("> data file list: {filtered_dir_content}".format(filtered_dir_content=str(filtered_dir_content)))
    verbosePrint("")
    POOL_alldata_dict = {}
    POOL_IS_dict = {}
    for path in filtered_dir_content:
        filename = str(os.path.basename(path))
        barcode = ".".join((filename.split("."))[:2])
        verbosePrint("> Processing {filename}, barcode={barcode} ...".format(filename=str(filename), barcode=str(barcode)))
        data_file_nested_list = loadFile(path, data_files_delimiter)
        alldata_dict, IS_dict = arrangeData(data_file_nested_list)
        POOL_alldata_dict[barcode] = alldata_dict
        POOL_IS_dict[barcode] = IS_dict
        # If verbose, maybe some details about file content should be printed (new function deepVerbosePrint?)
    verbosePrint(">>> Data Loaded!")
    return POOL_alldata_dict, POOL_IS_dict


#++++++++++++++++++++++++++++++++++++++ MAIN and TEST ++++++++++++++++++++++++++++++++++++++#

if __name__ == "__main__":
    
#    # Test vars - Gemini
#    ground_dir = "/opt/NGS/results"
#    DISEASE = "AssayValidation"
#    PATIENT = "CEMJY"
#    POOL = "LANE_1"
#    data_files_delimiter = '\t'
#    data_files_name_filter = ".randomBC.tsv"  # !!!
#
#    # Test code - Gemini
#    # loadDataFiles calls buildInputPath and then loops over files returned by listDir
#    # Each file is loaded by loadFile and arranged as nested dictionary by arrangeData: alldata_dict, IS_dict
#    # Each file is associated to a barcode, so results of the whole POOL are aggregated in POOL_alldata_dict, POOL_IS_dict, increasing nesting with 'barcode' keys.
#    POOL_alldata_dict, POOL_IS_dict = loadDataFiles(ground_dir, DISEASE, PATIENT, POOL, data_files_name_filter, data_files_delimiter)
    
    # Test vars - Local
    delimiter = '\t'
    path = "/home/stefano/Desktop/RandomBC_matrix_development/test_input/data"
    data_files_name_filter = ".randomBC.tsv"  # !!!
    filtered_dir_content = listDir(path, name_filter=data_files_name_filter)
    # Test code - Local
    verbosePrint("\n>>> Loading data ...")
    verbosePrint("> path: {path}".format(path=str(path)))
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
        data_file_nested_list = loadFile (path, delimiter)
        alldata_dict, IS_dict = arrangeData(data_file_nested_list)
        POOL_IS_dict[barcode] = IS_dict
        POOL_alldata_dict[barcode] = alldata_dict
    verbosePrint(">>> Data Loaded!")








 