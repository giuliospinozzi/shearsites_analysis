#!/usr/bin/python

import argparse, os, sys

header = """

+--------------------------------------------------+
 Author:    Stefano Brasca
 Date:    December 2015
 Contact:    brasca.stefano@hsr.it
 Version:    1.0
+--------------------------------------------------+

 Description
  - This is a pipeline plugin for Shear Site Random Barcode data build.
  
 Note:
  - the script needs files like *.shearsites.tsv (usually
    *.shearsites.ISfixed.LENGTHfixed.tsv) and *.randomBC.fasta files.

 Steps
    1. Load .shearsites.tsv (or as supplied through --origFile_nameEnd)
       from cwd (or as supplied by --folder)
    2. Load .randomBC.fasta (or as supplied through --BCfasta_nameEnd)
       from cwd (or as supplied by --folder)
    3. Couple shearsite data with random barcode sequences, through headers
    4. Write results in *.shearsites.ISfixed.LENGTHfixed.randomBC.tsv files
       (or as suppliedby --outFile_nameEnd)
""" 

description = """
This is a pipeline plugin for Shear Site Random Barcode data build.
"""

usage_example = ""

# print header
parser = argparse.ArgumentParser(usage = usage_example, epilog = "[ hSR-TIGET - Vector Integration Core - Bioinformatics ] \n", description = description)
parser.add_argument('--folder', dest="folder", help="path of the folder to operate into", action="store", default=os.getcwd(), required=False)
parser.add_argument('--origFile_nameEnd', dest="origFile_nameEnd", help="string composed by suffixes + extension, to identify original files", action="store", default=".shearsites.ISfixed.LENGTHfixed.tsv", required=False)
parser.add_argument('--BCfasta_nameEnd', dest="BCfasta_nameEnd", help="string composed by suffixes + extension to identify fasta files containing random BarCodes", action="store", default=".randomBC.fasta", required=False)
parser.add_argument('--outFile_nameEnd', dest="outFile_nameEnd", help="string to name out files, used as suffixes + extension", action="store", default=".shearsites.ISfixed.LENGTHfixed.randomBC.tsv", required=False)
#parser.add_argument('--an_arg', dest="an_arg", help="some help", action="store", required=True)
#parser.add_argument('--another_arg', dest="another_arg", help="some help again", action="store", default="default_value", required=False)
args = parser.parse_args()



#########################################################################################
### FUNCTIONS
#########################################################################################

#def checkArgs(args):
#    """
#    Check to be developed
#    """
#    if False:
#        print "\n[AP]\tAn error occurred\n\tExit\n"
#        sys.exit()

def pathList(folder_path, name_filter=None):
    """
    Returns the content (abs path list) of input folder, filtered
    by name_filter substring, if not None.
    """
    dir_content = []
    for filename in os.listdir(os.path.normpath(folder_path)):
        dir_content.append(os.path.abspath(filename))
    if name_filter is None:
        return dir_content
    else:
        filtered_dir_content = []
        for line in dir_content:
            if name_filter in line:
                filtered_dir_content.append(line)
        return filtered_dir_content
    
    
def loadFile(path, split_data='\t', end_of_line='\n'):
    """
    Load file as nested list
    """
    row_list = []
    with open(path) as in_file:
        for row in in_file:
            split_row = row.split(split_data)
            split_row[-1] = split_row[-1].replace(end_of_line, "")
            row_list.append(split_row)
        # Test Print: loading visual check
        #print "\nData loaded from: ", path
        #print "Total N_row =", len(row_list)
        #if len(row_list) != 0:
            #print "Total N_col =", len(row_list[0])
            #print "Visual check:"
            #print row_list[0][0], row_list[0][1], "[...]", row_list[0][-1]
            #if len(row_list) > 2:
                #print "[...]"
            #if len(row_list) > 1:
                #print row_list[-1][0], row_list[-1][1], "[...]", row_list[-1][-1]
    return row_list


def loadFASTAFile(path):
    """
    Load fasta file as nested list: [[header0, seq0], [header1, seq1], ...]
    """
    with open(path, 'r') as in_file:
        row_list = []
        source_file_lines = in_file.readlines()
        n_of_seq = len(source_file_lines)
        if n_of_seq == 0:
            return row_list
        elif n_of_seq%2 != 0:
            print "\n[ERROR] {path} FASTA file has an odd number of lines!\n".format(path=str(path))
        # Loop over source file lines
        header = None
        for source_file_line in source_file_lines:
            if '>' in source_file_line:
                tmp_header = source_file_line.rstrip()[1:]
                header = tmp_header.split(" ")[0]   # Test in Spyder
            else:
                if header:
                    sequence = source_file_line.rstrip()
                    row_list.append([header, sequence])
                    header = None
                else:
                    print "\n[ERROR] {path} FASTA file is weird!\n".format(path=str(path))
    return row_list
   

def buildData(original_file, BCfasta_file, out_file, out_file_separator='\t', out_file_end_of_line='\n'):
    """
    Couple random barcode data in BCfasta_file with data in original_file, through headers.
    Results in out_file. TSV format like input ones, with last column appended.
    """
    file_dict = {}
    fasta_dict = {}
    headers = []
    for line_as_list in loadFile(original_file):
        header = line_as_list[0]
        headers.append(header)
        file_dict[header] = line_as_list   # the whole line
    for line_as_list in loadFASTAFile(BCfasta_file):
        header = line_as_list[0]
        seq = line_as_list[1]
        fasta_dict[header] = seq   # the seq
    # check: headers supposed to be not duplicated
    if (set(file_dict.keys()) != set(fasta_dict.keys())):
        print "\n[ERROR] cannot build data with {original_file} and {BCfasta_file}!\n".format(original_file=str(original_file), BCfasta_file=str(BCfasta_file))
        shearsite_file = set(file_dict.keys())
        fasta_file = set(fasta_dict.keys())
        s = shearsite_file.difference(fasta_file)
        f = fasta_file.difference(shearsite_file)
        print "        * headers in shearsite file not in fasta file: ", sorted(list(s))
        print "        * headers in fasta file not in shearsite file: ", sorted(list(f))
        return 0
    # Compute outFile_lines
    outFile_lines = []
    for header in headers:
        out_line = file_dict[header] + [fasta_dict[header]]
        outFile_lines.append(out_line)
    # Write out_file
    with open(out_file, 'w') as out_stream:
        for line in outFile_lines:
            out_stream.write(out_file_separator.join(line) + out_file_end_of_line)
    return 0



#########################################################################################
### MAIN
#########################################################################################

def main():
    """
    Main part of the program.
    """
    #checkArgs(args)
    
    # Select input files
    original_files_to_process = sorted(pathList(args.folder, args.origFile_nameEnd))
    BC_files_to_process = sorted(pathList(args.folder, args.BCfasta_nameEnd))
    # Prepare output paths
    out_file_absPaths = [a.replace(args.origFile_nameEnd, args.outFile_nameEnd) for a in original_files_to_process]
    
    # Test print
    #print original_files_to_process
    #print BC_files_to_process
    #print out_file_absPaths
    
    # Build final files
    for original_file, BC_file, out_file in zip(original_files_to_process, BC_files_to_process, out_file_absPaths):
        buildData(original_file, BC_file, out_file)
        
    return 0
    
    
# sentinel
if __name__ == "__main__":
    main()