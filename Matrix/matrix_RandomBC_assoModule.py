# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 13:57:28 2016

@author: stefano
"""


#++++++++++++++ Requested Package(s) Import +++++++++++++++#
import os, sys
import csv
import re


#++++++++++++++++++++++ Global Vars +++++++++++++++++++++++#
verbose = True


#+++++++++++++++++++++++++++++++++++++++ FUNCTIONS +++++++++++++++++++++++++++++++++++++++#

def verbosePrint(x, verbose=verbose):
    if verbose:
        print x

def humanSorted(l):
    def tryint(s):
        try:
            return int(s)
        except:
            return s
    def alphanum_key(s):
        return [ tryint(c) for c in re.split('([0-9]+)', s) ]
    def sort_nicely(l):
        return sorted(l, key=alphanum_key)
    return sort_nicely(l)
    

def generateFilePath (filename, folder=os.getcwd()):
    """
    Purpose: WRITE A FILE
    
    IN: 
    filename - ... (trying to crack with a path raise an error)
    folder - ... complete path (e.g. starting with '/') are checked and taken as they are.
                 words or partial path ('foldername' or 'foldername/subfolder') will be
                 taken as cwd/...
    OUT:
    complete_path - complete path from folder+filename (checked)
    """
    # Check folder and create if it doesn't exist
    try:
        folder = os.path.abspath(os.path.normpath(folder))   ### good mod, backward compatible!
    except Exception, err_message:
        print "\n[ERROR] '{folder}' is not formatted as valid output path!".format(folder=str(folder))
        print "os.path.normpath returned: ", err_message
        sys.exit("\n[QUIT]\n")
    if os.path.isfile(folder):
        print "\n[ERROR] Output path must be a folder path, not a file! Your input: '{folder}'".format(folder=str(folder))
        sys.exit("\n[QUIT]\n")
    if not os.path.exists(folder):
        try:
            os.makedirs(folder)
        except Exception, err_message:
            print "\n[ERROR] '{folder}' is not a valid folder path or you don't have sufficient privileges to create it!".format(folder=str(folder))
            print "os.makedirs returned: ", err_message
            sys.exit("\n[QUIT]\n")
    if not os.access(folder, os.W_OK):
        print "\n[ERROR] Can't write anything in {folder}'".format(folder=str(folder))
        sys.exit("\n[QUIT]\n")
    # Cast filename
    filename = str(filename)
    # Create complete path and warn if exists!
    complete_path = os.path.join(folder, filename)
    if os.path.isfile(complete_path):
        print "\n[WARNING] {complete_path} already exists!".format(complete_path=str(complete_path))
    # Simulate writing
    try:
        filehandle = open(complete_path, 'w')
        filehandle.close()
    except IOError:
        print "\n[ERROR] '{complete_path}' is not a valid file path or you don't have sufficient privileges to write there!".format(complete_path=str(complete_path))
        sys.exit("\n[QUIT]\n")
    return complete_path
    
    
def getFilePath (filename, folder=os.getcwd()):
    """
    Purpose: READ A FILE
    
    IN: 
    filename - ... (trying to crack with a path raise an error)
    folder - ... complete path (e.g. starting with '/') are checked and taken as they are.
                 words or partial path ('foldername' or 'foldername/subfolder') will be
                 taken as cwd/...
    OUT:
    complete_path - complete path from folder+filename (checked)
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
        print "\n[ERROR] You have not read privileges in folder='{folder}'".format(folder=str(folder))
        sys.exit("\n[QUIT]\n")
    if not os.access(folder, os.W_OK):
        print "\n[WARNING] You have not write privileges in folder='{folder}'".format(folder=str(folder))
    # Cast filename
    filename = str(filename)
    # Create complete path and check if exists!
    complete_file_path = os.path.join(folder, filename)
    if not os.path.isfile(complete_file_path):
        print "\n[ERROR] complete_file_path='{complete_file_path}' does not exist!".format(complete_file_path=str(complete_file_path))
        sys.exit("\n[QUIT]\n")
    # Test the readability of complete_file_path
    if not os.access(complete_file_path, os.R_OK):
        print "\n[ERROR] You have not read privileges for complete_file_path='{complete_file_path}'".format(complete_file_path=str(complete_file_path))
        sys.exit("\n[QUIT]\n")
    return complete_file_path


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
    
    
def arrangeAssoData(asso_file_nested_list):
    """
    Purpose: ARRANGE NESTED LIST AS NESTED DICT
    IN:
    asso_file_nested_list
    OUT:
    asso_dict = {'masterkey0': {linedict0}, 'masterkey1': {linedict1}, ...}
                linedict = {'sub-key0': cell_of_field0, 'sub-key1': cell_of_field1, ...}
    Notes:
    master-keys are generated by generateMasterKey func.
    sub-dict and sub-keys are generated by generateLineDict.
    """
    def generateMasterKey(line, master_key_indexes=[0], join_char='_'):
        """
        master_key_indexes = [0] means master_key -> First field
        of assoFile (barcode!)
        """
        master_key_items = []
        master_key_items_append = master_key_items.append
        for i in master_key_indexes:
            master_key_items_append(line[i])
        return join_char.join(master_key_items)
    def generateLineDict(line):
        """
        from line (list) to line_dict (dict with indexes as keys)
        """
        line_dict = {}
        for i in range(0,len(line)):
            line_dict[i] = line[i]
        return line_dict
    
    asso_dict = {}
    for line in asso_file_nested_list:
        master_key = generateMasterKey(line)
        line_dict = generateLineDict(line)
        asso_dict[master_key] = line_dict
    if len(asso_file_nested_list) != len(asso_dict):
        print "\n[ERROR] arrangeAssoData generated not unique master-keys for asso_dict! Check Association File contents and compare to generateMasterKey sub-function!"
        sys.exit("\n[QUIT]\n")
    return asso_dict
    

def loadAssoFile(asso_file_name, asso_folder, asso_delimiter):
    # Load asso_file
    asso_file_complete_path = getFilePath (asso_file_name, folder=asso_folder)
    asso_file_nested_list = loadFile (asso_file_complete_path, delimiter=asso_delimiter)
    # Check asso_file (must be not empty) and fields (must have a 'square' structure)
    first_line = None
    n_fields = None
    fields_consistence = True
    try:
        first_line = asso_file_nested_list[0]
        n_fields = len(first_line)
    except:
        print "\n[ERROR] asso_file_complete_path='{asso_file_complete_path}' is empty!".format(asso_file_complete_path=str(asso_file_complete_path))
        sys.exit("\n[QUIT]\n")
    if not (n_fields > 0):
        print "\n[ERROR] asso_file_complete_path='{asso_file_complete_path}' has no columns!".format(asso_file_complete_path=str(asso_file_complete_path))
        sys.exit("\n[QUIT]\n")
    for line in asso_file_nested_list:
        if len(line) != n_fields:
            fields_consistence = False
            n_fields = max([len(line), n_fields])
    if fields_consistence is False:
        ### TO DO: Here implement a function to let asso_file_nested_list be a "square" filled with "something" where needed!
        ### Remembrer to update the final summary!
        print "\n[ERROR] asso_file_complete_path='{asso_file_complete_path}' has not a regular structure (is not a 'square')!".format(asso_file_complete_path=str(asso_file_complete_path))
        sys.exit("\n[QUIT]\n")
    # Arrange data
    asso_dict = arrangeAssoData(asso_file_nested_list)
    # Print a Summary if verbose
    verbosePrint("\n>>> Association File loaded!")
    verbosePrint("> path = {path}".format(path=str(asso_file_complete_path)))
    master_keys = humanSorted(asso_dict.keys())
    verbosePrint("> n lines = {n_lines}".format(n_lines=str(len(master_keys))))
    sub_keys = humanSorted([str(x) for x in asso_dict[master_keys[0]].keys()])
    verbosePrint("> n fields = {n_fields}".format(n_fields=str(len(sub_keys))))
    verbosePrint("> master-keys:")
    verbosePrint(master_keys)
    verbosePrint("> sub-keys:")
    verbosePrint(sub_keys)
    verbosePrint("")
    # Return
    return asso_dict


#++++++++++++++++++++++++++++++++++++++ MAIN and TEST ++++++++++++++++++++++++++++++++++++++#

if __name__ == "__main__":
    # Input vars
    asso_folder = "/home/stefano/Desktop/RandomBC_matrix_development/test_input/asso"  # /opt/applications/scripts/isatk/elements/association
    asso_file_name = "asso.assayvalidation.lane1.tsv"
    asso_delimiter = '\t'
    # loadAssoFile calls getFilePath, then loadFile.
    # Loaded data are arranged by arrangeAssoData, which defines asso_dict structure through generateMasterKey and generateLineDict
    asso_dict = loadAssoFile(asso_file_name, asso_folder, asso_delimiter)











