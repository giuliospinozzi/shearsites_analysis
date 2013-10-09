#!/usr/bin/python

import argparse, os

header = """

+--------------------------------------------------+
 Author:	Giulio Spinozzi
 Date:		September 2013
 Contact:	spinozzi.giulio@hsr.it
 Version:	2.0
+--------------------------------------------------+

 Description
  - This is a pipeline plugin for the Shear Sites identification.
  
 Note:
  - the script needs some files: the initial R2 BED (LTR96.LC41.r2.bed), the final R1 BED (LTR96.LC41.sorted.rel.pg.dupflag.tagfilter.cigarmd.bed)

 Steps
	1. loads the input files
	2. creates two dictionary with ID the header of R1 and R2
	3. makes the comparison between them
	4. create the list of the shear sites (candidate or real)
""" 

description = "This is a pipeline plugin for the Shear Sites identification."

usage_example = """
Examples of usage:
APP --bed1 LTR96.LC41.sorted.rel.pg.dupflag.tagfilter.cigarmd.bed --bed2 LTR96.LC41.r2.bed
"""

# print header

parser = argparse.ArgumentParser(usage = usage_example, epilog = "[ hSR-TIGET - Vector Integration Core - Bioinformatics ] \n", description = description)
parser.add_argument('--bed1', dest="bedfile_r1", help="BED file of R1 to process. No default option.", action="store", required=True)
parser.add_argument('--bed2', dest="bedfile_r2", help="BED file of R2 to process. No default option.", action="store", required=True)
parser.add_argument('--all', dest="all", help="With this option the output file contains all the candidate shear sites", action="store_true", default=False, required=False)
args = parser.parse_args()



#########################################################################################
### FUNCTIONS
#########################################################################################

def checkArgs(args):
	"""
	Check file path
	"""
	if not os.path.isfile(args.bedfile_r1) or not os.path.isfile(args.bedfile_r2):
		print "\n[AP]\tError while reading files: no valid paths.\n\tExit\n"
		sys.exit()


# Check if R2 bed has the /2!!
def compare_R1R2_BEDs(bed_r1, bed_r2):
	with open(bed_r1) as input:
		a = zip(*(line.strip().split('\t') for line in input))
		set_r1 = set()
		for k in a[3]:
			b = k.split('/1')
			b = b[0]
			set_r1.add(b)
	print "\n[AP]\tR1 set loaded!"
	with open(bed_r2) as input:
		c = zip(*(line.strip().split('\t') for line in input))
		set_r2 = set()
		for j in c[3]:
			b = j.split('/2')
			b = b[0]
			set_r2.add(b)
	print "[AP]\tR2 set loaded!"
	count = 0
	for x in set_r1:
		if x in set_r2:
			count = count + 1
	print "\n[AP]\tCandidate Shear Sites: " + str(count)


# Generate the two dictionary, for R1 and R2
def generate_BEDs_dict(bed_r1, bed_r2):
	f_r1 = open(bed_r1,"r")
	f_r2 = open(bed_r2,"r")
	dict_r1 = {}
	dict_r2 = {}
	array = []
	for line in f_r1:
		array.append(line)
		string_splitted = line.split('\t')
		string_splitted[5] = string_splitted[5].split('\n')[0]
		key = string_splitted[3].split('/1')
		key = key[0]
		dict_r1[key] = string_splitted[0], string_splitted[1], string_splitted[2], string_splitted[4], string_splitted[5]
	# print "\nR1 Keys: " 
	# print dict_r1.keys()
	for line in f_r2:
		array.append(line)
		string_splitted = line.split('\t')
		string_splitted[5] = string_splitted[5].split('\n')[0]
		key = string_splitted[3].split('/2')
		key = key[0]
		quality = int(string_splitted[4])
		if dict_r2.has_key(key):
			if int(dict_r2[key][3]) < quality:
				dict_r2[key] = string_splitted[0], string_splitted[1], string_splitted[2], string_splitted[4], string_splitted[5]
		else:
			dict_r2[key] = string_splitted[0], string_splitted[1], string_splitted[2], string_splitted[4], string_splitted[5]
	# print "\nR2 Keys: "
	# print dict_r2.keys()

	return dict_r1, dict_r2


# Shear site identification, it generates the output file (all the candidate shear sites or only the real ones)
def shearSites_identification(dict_r1, dict_r2, outfilename, parameter):
	# write python dict to a file

	if parameter == True:
	# All dataset, data not filtered (a bed file with all reads, R1-R2 as they are)
		f_out = open(outfilename, 'w')
		for key in dict_r1.keys():
			values_r1 = dict_r1[key]
			if dict_r2.has_key(key):
				values_r2 = dict_r2[key]
				f_out.write(key + '\t' + values_r1[0] + '\t' + values_r1[1] + '\t' + values_r1[2] + '\t' + values_r1[3] + '\t' + values_r1[4] + '\t' + values_r2[0] + '\t' + values_r2[1] + '\t' + values_r2[2] + '\t' + values_r2[3] + '\t' + values_r2[4] + '\n')
			else:
				f_out.write(key + '\t' + values_r1[0] + '\t' + values_r1[1] + '\t' + values_r1[2] + '\t' + values_r1[3] + '\t' + values_r1[4] + '\n')
	else:
	# Real shear sites
		f_out = open(outfilename, 'w')
		for key in dict_r1.keys():
			values_r1 = dict_r1[key]
			if dict_r2.has_key(key):
				values_r2 = dict_r2[key]
				len_plus = abs(int(values_r1[1]) - int(values_r2[2]))
				len_minus = abs(int(values_r1[2]) - int(values_r2[1]))
				if values_r1[0] == values_r2[0]:
					if (values_r1[4] == '+' and values_r2[4] == '-' and len_plus < 1200):
						f_out.write(key + '\t' + values_r1[0] + '\t' + values_r1[1] + '\t' + values_r1[2] + '\t' + values_r1[3] + '\t' + values_r1[4] + '\t' + values_r2[0] + '\t' + values_r2[1] + '\t' + values_r2[2] + '\t' + values_r2[3] + '\t' + values_r2[4] + '\t' + str(len_plus) + '\n')
					elif (values_r1[4] == '-' and values_r2[4] == '+') and (len_minus < 1200):
						f_out.write(key + '\t' + values_r1[0] + '\t' + values_r1[1] + '\t' + values_r1[2] + '\t' + values_r1[3] + '\t' + values_r1[4] + '\t' + values_r2[0] + '\t' + values_r2[1] + '\t' + values_r2[2] + '\t' + values_r2[3] + '\t' + values_r2[4] + '\t' + str(len_minus) + '\n')

	f_out.close()


# # This method generates the output file format for R analysis with Berry Model (sonicLength)
# def generateOutputForR(path):
# 	# The system call produces a txt file with chr, locus of integration (start), strand and length of fragment for each tsv.
# 	# Unique tuples and non unique
# 	call = "for k in $(ls *.shearsites.tsv); do awk '{print $2\" \"$3\" \"$6\"\t\"$12}' %(path)s$k | sort > %(path)s${k:0:21}.txt ;done" %{	'path': path,
# 	}
# 	call_uniq = "for k in $(ls *.shearsites.tsv); do awk '{print $2\" \"$3\" \"$6\"\t\"$12}' %(path)s$k | sort | uniq > %(path)s${k:0:21}.uniq.txt ;done" %{	'path': path,
# 	}
# 	os.system(call)
# 	os.system(call_uniq)



#########################################################################################
### MAIN
#########################################################################################

def main():
	"""
	Main part of the program.
	"""
	#first check args and file paths
	checkArgs(args)
	print "[AP]\tChecked inputs, now acquiring data"

	bed_r1 = args.bedfile_r1
	bed_r2 = args.bedfile_r2

	if os.path.getsize(bed_r1) > 0 and os.path.getsize(bed_r2) > 0:
		f_name = bed_r1.split(".")
		outfilename = f_name[0] + '.' + f_name[1] + '.shearsites.tsv'

		compare_R1R2_BEDs(bed_r1, bed_r2)

		dict_r1, dict_r2 = generate_BEDs_dict(bed_r1, bed_r2)

		parameter = args.all
		shearSites_identification(dict_r1, dict_r2, outfilename, parameter)

		# generateOutputForR("/home/giulio/Dropbox/tmp/")
	else:
		print "\n[AP]\tEmpty input file: no output!!"

	print "\n[AP]\tTask Finished, closing.\n"


# sentinel
if __name__ == "__main__":
    main()