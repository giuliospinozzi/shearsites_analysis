#!/usr/bin/python

import argparse, os, sys

header = """

+--------------------------------------------------+
 Author:    Stefano Brasca
 Date:    August 2015
 Contact:    brasca.stefano@hsr.it
 Version:    1.0
+--------------------------------------------------+

 Description
  - This is a pipeline plugin for the Shear Site aggregation around ISs.
  
 Note:
  - the script needs *.shearsites.raw.sort.uniq.bed, ISrange.bed,
    IScoordinate.bed

 Steps
    1. Load ISrange.bed (or as supplied through --range_filename)
       from cwd (or as supplied by --folder)
    1. Load IScoordinate.bed (or as supplied through --IS_filename)
       from cwd (or as supplied by --folder)
    2. Load *.shearsites.raw.sort.uniq.bed files (or as supplied through --infiles_nameEnd)
       from cwd (or as supplied by --folder)
    3. Fix mapping coordinates: entries mapping inside ranges provided by ISrange.bed
       are shifted to the related IS positions provided by IScoordinate.bed
    4. Write results in *.shearsites.ISfixed.bed files (or as supplied
       by --outFile_nameEnd)
""" 

description = """
This is a pipeline plugin for the Shear Site aggregation around ISs.
"""

usage_example = ""

# print header
parser = argparse.ArgumentParser(usage = usage_example, epilog = "[ hSR-TIGET - Vector Integration Core - Bioinformatics ] \n", description = description)
parser.add_argument('--folder', dest="folder", help="path of the folder to operate into", action="store", default=os.getcwd(), required=False)
parser.add_argument('--range_filename', dest="range_filename", help="alternative name for ISrange.bed input file", action="store", default="ISrange.bed", required=False)
parser.add_argument('--IS_filename', dest="IS_filename", help="alternative name for IScoordinate.bed input file", action="store", default="IScoordinate.bed", required=False)
parser.add_argument('--inFiles_nameEnd', dest="inFiles_nameEnd", help="string composed by suffixes + extension, to identify input bed files", action="store", default=".shearsites.raw.sort.uniq.bed", required=False)
parser.add_argument('--outFile_nameEnd', dest="outFile_nameEnd", help="string to name out files, used as suffixes + extension", action="store", default=".shearsites.ISfixed.bed", required=False)
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
    

def build_ISdict(IS_file, range_file):
    """
    """
    IS_dict = {}
    for IS_line, range_line in zip(loadFile(IS_file), loadFile(range_file)):
        chromosome = range_line[0]
        startLocus = range_line[1]
        endLocus = range_line[2]
        strand = range_line[-1]
        integrationLocus = IS_line[1]
        for n in range(int(startLocus), int(endLocus)+1):
            key = "_".join([chromosome, str(n), strand])
            IS_dict[key] = integrationLocus
        # soft check # comment after test!
        if ((chromosome != IS_line[0]) or (strand != IS_line[-1])):
            print "\n[ERROR]\t{IS_file} and {range_file} are weird or not properly paired! Chr and/or strand don't match!\n".format(IS_file=str(IS_file), range_file=str(range_file))
            return False
        if ((int(integrationLocus) > int(endLocus)) or (int(integrationLocus) < int(startLocus))):
            print "\n[ERROR]\t{IS_file} and {range_file} are weird or not properly paired! startLocus <= integrationLocus <= endLocus condition is violated\n".format(IS_file=str(IS_file), range_file=str(range_file))
            return False
    return IS_dict
        

def applyCorrection(in_file, out_file, IS_dict, out_file_separator='\t', out_file_end_of_line='\n'):
    """
    
    """
    # Compute outFile_lines
    outFile_lines = []
    i = 0
    for in_line in loadFile(in_file):
        chromosome = in_line[0]
        locus = in_line[1]
        strand = in_line[-2]
        fragLength = in_line[-1]
        key = "_".join([chromosome, locus, strand])
        out_line = [chromosome, IS_dict[key], IS_dict[key], "", "", strand, fragLength]
        outFile_lines.append(out_line)
        # soft check
        i+=1
        if locus != in_line[2]:
            print "\n[ERROR]\tline {i} in {in_file} has different start and end locus!\n".format(in_file=str(in_file), i=str(i))
            print "\tprint line as list: ", in_line
            print "\tSKIP THIS FILE!\n"
            return 0
        
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
    ISfile_path = pathList(args.folder, args.IS_filename)[0]
    rangeFile_path = pathList(args.folder, args.range_filename)[0]
    in_file_absPaths = sorted(pathList(args.folder, args.inFiles_nameEnd))
    # Prepare output paths
    out_file_absPaths = [a.replace(args.inFiles_nameEnd, args.outFile_nameEnd) for a in in_file_absPaths]
    
    # Test print
    #print ISfile_path, " == ", ISfile_path[0]
    #print rangeFile_path, " == ", rangeFile_path[0]
    #print in_file_absPaths
    #print out_file_absPaths
    
    # Load IS data and Range data in a unique struncture
    IS_dict = build_ISdict(ISfile_path, rangeFile_path)
    
    # Create corrected files
    if IS_dict:
        for in_file, out_file in zip(in_file_absPaths, out_file_absPaths):
            applyCorrection(in_file, out_file, IS_dict)
    else:
        print "Unable to build range and IS data for aggregation\n\t[EXIT]\n"
        sys.exit()
        
    return 0
    
    
# sentinel
if __name__ == "__main__":
    main()