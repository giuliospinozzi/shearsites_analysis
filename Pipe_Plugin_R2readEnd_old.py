#!/usr/bin/python

import argparse, os

header = """

+--------------------------------------------------+
 Author:	Giulio Spinozzi
 Date:		August 2013
 Contact:	spinozzi.giulio@hsr.it
 Version:	1.0
+--------------------------------------------------+

 Description
  - This is a pipeline plugin for update the read score start-end in the final BED file (after cigar filter).
  
 Note:
  - the script needs some files: the initial R2 BED (LTR96.LC41.r2.bed), the final R1 BED (LTR96.LC41.sorted.rel.pg.dupflag.tagfilter.cigarmd.bed)
  - All R1 entries in the final BED must be included in the R2 initial BED

 Steps
	1. loads the input files
	2. creates two dictionary with ID the header of R1 and R2
	3. makes the comparison between them
	4. updates the cigar BED file with the R2 start-end locus
""" 

description = "This is a pipeline plugin for update the read score start-end in the final BED file (after cigar filter)."

usage_example = """
Examples of usage:
APP --bed1 LTR96.LC41.sorted.rel.pg.dupflag.tagfilter.cigarmd.bed --bed2 LTR96.LC41.r2.bed -o LTR96.LC41.sorted.rel.pg.dupflag.tagfilter.cigarmd.sort.product.bed --stdOut true
"""

print header

parser = argparse.ArgumentParser(usage = usage_example, epilog = "[ hSR-TIGET - Vector Integration Core - Bioinformatics ] \n", description = description)
parser.add_argument('--bed1', dest="bedfile_r1", help="BED file of R1 to process. No default option.", action="store", required=True)
parser.add_argument('--bed2', dest="bedfile_r2", help="BED file of R2 to process. No default option.", action="store", required=True)
parser.add_argument('-o', '--outfilename', dest="outfilename", help="The output BED file. No default option.", action="store", required=True)
parser.add_argument('--stdOut', dest="stdOut", help="Flag for use the standard output for generate the final BED. No default option.", action="store_true", default=False)
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


def compare_R1R2_BEDs(bed_r1, bed_r2):
	with open(bed_r1) as input:
		a = zip(*(line.strip().split('\t') for line in input))
		set_r1 = set()
		for k in a[3]:
			b = k.split('/1')
			b = b[0]
			set_r1.add(b)
	print "R1 set loaded!"
	with open(bed_r2) as input:
		c = zip(*(line.strip().split('\t') for line in input))
		set_r2 = set()
		for j in c[3]:
			set_r2.add(j)
	print "R2 set loaded!"
	print "Starting comparison..."
	count = 0
	for x in set_r1:
		if x in set_r2:
			count = count + 1
	print count
	print "Comparison Done!"


def generate_BED(bed_r1, bed_r2, outfilename):
	f_r1 = open(bed_r1,"r")
	f_r2 = open(bed_r2,"r")
	dict_r1 = {}
	dict_r2 = {}
	array = []
	for line in f_r1:
	    array.append(line)
	    str = line.split('\t')
	    str[5] = str[5].split('\n')[0]
	    dict_r1[str[3]] = str[0], str[1], str[2], str[4], str[5]
	for line in f_r2:
	    array.append(line)
	    str = line.split('\t')
	    str[5] = str[5].split('\n')[0]
	    dict_r2[str[3]+'/1'] = str[0], str[1], str[2], str[4], str[5]

	# write python dict to a file
	f_out = open(outfilename, 'w')
	for key, value in dict_r1.iteritems():
		values_r2 = dict_r2[key]
		if values_r2[3]=='0':
			f_out.write(value[0] + '\t' + value[1] + '\t' + value[2] + '\t'+ key + '\t' + value[3] + '\t' + value[4] + '\t' + '0' + '\n')
		else:
			if (value[4]=='+' and values_r2[4]=='-'):
				f_out.write(value[0] + '\t' + value[1] + '\t' + values_r2[2] + '\t'+ key + '\t' + value[3] + '\t' + value[4] + '\t' + '1' + '\n')
			if (value[4]=='-' and values_r2[4]=='+'):
				f_out.write(value[0] + '\t' + values_r2[1] + '\t' + value[2] + '\t'+ key + '\t' + value[3] + '\t' + value[4] + '\t' + '1' + '\n')

	f_r1.close()
	f_r2.close()
	f_out.close()



#########################################################################################
### MAIN
#########################################################################################

def main():
	"""
	Main part of the program.
	"""
	#first check args and file paths
	checkArgs(args)
	print "[AP]\tChecked inputs, now acquiring data, BED, FASTQ without LTR and FASTQ with LTR"

	bed_r1 = args.bedfile_r1
	bed_r2 = args.bedfile_r2
	outfilename = args.outfilename
	outfilename2 = outfilename[:-3] + 'final.bed'

	standardOutput = args.stdOut

	compare_R1R2_BEDs(bed_r1, bed_r2)

	generate_BED(bed_r1, bed_r2, outfilename)

	if standardOutput=='true':
		call = "cat %(outfilename)s | bedtools sort > %(outfilename2)s" 	%{	'outfilename': outfilename,
			'outfilename2': outfilename2,
		 }
		os.system(call)
	else:
		call = "bedtools sort -i %(outfilename)s > %(outfilename2)s" 	%{	'outfilename': outfilename,
			'outfilename2': outfilename2,
		 }
		os.system(call)

	print "\n[AP]\tTask Finished, closing.\n"


# sentinel
if __name__ == "__main__":
    main()