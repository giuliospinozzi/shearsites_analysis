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
  - This is a pipeline plugin for the Shear Site length correction, after
    aggregation around ISs.
  
 Note:
  - the script needs *.shearsites.raw.sort.uniq.bed and 
    *.shearsites.ISfixed.bed files. Such files MUST BE PAIRED LINE BY LINE.
    Controls NOT YET IMPLEMENTED, only 'soft checks'.

 Steps
    1. Load *.shearsites.raw.sort.uniq.bed (or as supplied through --origFile_nameEnd)
       from cwd (or as supplied by --folder)
    2. Load .shearsites.ISfixed.bed (or as supplied through --ISfixedFile_nameEnd)
       from cwd (or as supplied by --folder)
    3. Fix Shear Site lengths by line-by-line comparison
    4. Write results in *.shearsites.ISfixed.LENGTHfixed.bed files (or as supplied
       by --outFile_nameEnd)
""" 

description = """
This is a pipeline plugin for the Shear Site length correction, after aggregation around ISs.
"""

usage_example = ""

# print header
parser = argparse.ArgumentParser(usage = usage_example, epilog = "[ hSR-TIGET - Vector Integration Core - Bioinformatics ] \n", description = description)
parser.add_argument('--folder', dest="folder", help="path of the folder to operate into", action="store", default=os.getcwd(), required=False)
parser.add_argument('--origFile_nameEnd', dest="origFile_nameEnd", help="string composed by suffixes + extension, to identify original files", action="store", default=".shearsites.raw.sort.uniq.bed", required=False)
parser.add_argument('--ISfixedFile_nameEnd', dest="ISfixedFile_nameEnd", help="string composed by suffixes + extension to identify IS fixed files", action="store", default=".shearsites.ISfixed.bed", required=False)
parser.add_argument('--outFile_nameEnd', dest="outFile_nameEnd", help="string to name out files, used as suffixes + extension", action="store", default=".shearsites.ISfixed.LENGTHfixed.bed", required=False)
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
    
    
def applyCorrection(original_file, ISfixed_file, out_file, out_file_separator='\t', out_file_end_of_line='\n'):
    """
    Fix Shear Site lengths by line-by-line comparison of original_file and ISfixed_file.
    Results in out_file. Bed format like input ones.
    """
    # Compute outFile_lines
    outFile_lines = []
    for original_line, IS_fixed_line in zip(loadFile(original_file), loadFile(ISfixed_file)):
        # soft check: only coordinates MAY be different
        if ((original_line[0] != IS_fixed_line[0]) or (original_line[-2] != IS_fixed_line[-2]) or (original_line[-1] != IS_fixed_line[-1])):
            print "\n[ERROR]\t{original_file} and {ISfixed_file} are properly paired!\n\tExit\n".format(original_file=str(original_file), ISfixed_file=str(ISfixed_file))
            sys.exit()
        # soft check: start-end coordinates must be equal within a file
        if ((original_line[1] != original_line[2]) or (IS_fixed_line[1] != IS_fixed_line[2])):
            print "\n[ERROR]\tstart-end coordinates don't match while processing {original_file} / {ISfixed_file}\n\tExit\n".format(original_file=str(original_file), ISfixed_file=str(ISfixed_file))
        # length correction
        if original_line[1] != IS_fixed_line[1]:
            diff = int(original_line[1]) - int(IS_fixed_line[1])
            fragLength = int(original_line[-1])
            if original_line[1] == '+':
                fragLength = fragLength + diff
            else:
                fragLength = fragLength - diff
            outFile_lines.append(IS_fixed_line[:-1] + [str(fragLength)])
            #Test print: changed lines
            #print outFile_lines[-1]
        else:
            outFile_lines.append(original_line)
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
    ISfixed_files_to_process = sorted(pathList(args.folder, args.ISfixedFile_nameEnd))
    # Prepare output paths
    out_file_absPaths = [a.replace(args.origFile_nameEnd, args.outFile_nameEnd) for a in original_files_to_process]
    
    # Test print
    #print original_files_to_process
    #print ISfixed_files_to_process
    #print out_file_absPaths
    
    # Create corrected files
    for original_file, ISfixed_file, out_file in zip(original_files_to_process, ISfixed_files_to_process, out_file_absPaths):
        applyCorrection(original_file, ISfixed_file, out_file)
        
    return 0
    
    
# sentinel
if __name__ == "__main__":
    main()