#!/usr/bin/python

import argparse, os
import MySQLdb

# Required
import rpy2.robjects as robjects
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import math

# Import sonicLength package (Rpy2 2.1+)
from rpy2.robjects.packages import importr
sonicLength = importr("sonicLength")

from scipy.stats import gaussian_kde

# Uncomment for Rpy2 2.0
# from rpy2.robjects import r
# r.library("sonicLength")


header = """

+--------------------------------------------------+
 Author:	Stefano Brasca, Giulio Spinozzi
 Date:		October 2013
 Contact:	brasca.stefano, spinozzi.giulio @hsr.it
 Version:	2.0
+--------------------------------------------------+

 Description
  - This is a wrapper with rpy2 for simulate the Berry's model in R.
  
 Note:
  - The script needs the dataset file (LTR60.LC90.shearsites.uniq.txt)

 Steps
	1. Loads the input file
	2. Calls the Berry's Model in R
	3. Saves the output of the model in Python structures
	4. Gives the output files with Phi and Theta
""" 

description = "This is a wrapper with rpy2 for simulate the Berry's model in R."

usage_example = """
Examples of usage:
APP --dataset LTR60.LC90.shearsites.uniq.txt
"""

# print header

parser = argparse.ArgumentParser(usage = usage_example, epilog = "[ hSR-TIGET - Vector Integration Core - Bioinformatics ] \n", description = description)
parser.add_argument('--dataset', dest="dataset_file", help="Dataset file to process in R. No default option.", action="store", required=True)
args = parser.parse_args()



#########################################################################################
### FUNCTIONS
#########################################################################################

def checkArgs(args):
	"""
	Check file path
	"""
	if not os.path.isfile(args.dataset_file):
		print "\n[AP]\tError while reading files: no valid paths.\n\tExit\n"
		sys.exit()


def generateDataset(data):
	# Load dataset
	dataset_file = open(data, "r")
	dataset = dataset_file.readlines()
	dataset_file.close()
	locations_list = []
	length_list = []
	for row in dataset:
		row = row.rstrip('\n')
		row_split = row.split("\t")
		locations_list.append(row_split[0])
		length_list.append(row_split[1])
	return locations_list, length_list


def fragmentsLengthPlot(length_phi,freq_phi,length_list,nameFile):
	length_phi_numbers = []
	for num in length_phi:
		length_phi_numbers.append(float(num))
	plt.figure()
	# Phi plot
	plt.plot(length_phi_numbers, freq_phi, 'r', hold=True, label="extimated distribution")
	length_list_numbers = []
	for num in length_list:
		length_list_numbers.append(float(num))
	binning = math.ceil(len(set(length_list_numbers))/2)
	# Plot (length-frequency) of input data
	plt.hist(length_list_numbers, bins=binning, normed=True, facecolor='green', alpha=0.3, hold=True, label="real distribution - histogram")
	density = gaussian_kde(length_list_numbers)
	xs = np.linspace(0,max(length_list_numbers)+25,len(set(length_list_numbers))*10+250)
	# Gaussian kde plot
	plt.plot(xs,density(xs), hold=True, label="real distribution - gaussian kde")
	fileName = nameFile + ".fragmentsLengthPlot.pdf"
   	plt.legend()
	# Labels
	plt.xlabel('fragments length')
	plt.ylabel('probability')
	plt.title(nameFile + ' Fragments length data')
   	plt.savefig(fileName, format="pdf")
   	return length_phi_numbers


def printThetaInfo(estimations_theta,locations_theta,nameFile):
	f_out = open(nameFile+".theta.tsv","w")
	dict_theta = {}
	for l in range(0,len(locations_theta)):
		dict_theta[locations_theta[l]] = estimations_theta[l]
		f_out.write(str(locations_theta[l]) + '\t' + str(estimations_theta[l]) + '\n')


def querySeqCount(host,user,passwd,db,db_table,destfile,nameFile):
	# System Call
	query = "mysql -h %(host)s -u %(user)s --password=%(passwd)s %(db)s --skip-column-names -e \"SELECT chr, integration_locus, strand, count(*) AS sequence_count FROM %(db_table)s WHERE tag='%(nameFile)s' GROUP BY chr, integration_locus,  strand ORDER BY chr ASC , integration_locus ASC\" > %(destfile)s " %{
     'host': host,
     'user': user,
     'passwd': passwd,
     'db': db,
     'db_table': db_table,
     'nameFile': nameFile,
     'destfile': destfile,
    }
	os.system(query)
	sequence_count = {}
	f_in = open(destfile, "r")
	for line in f_in:
		string_splitted = line.split('\t')
		value = string_splitted[3].split('\n')
		sequence_count["chr" + string_splitted[0] + " " + string_splitted[1] + " " + string_splitted[2]] = value[0]
	return sequence_count


def box_plot (real_data_list, extimated_data_list, dataset_name):
	"""
	dataset_name - string: dataset ID
	real_data_list - list of sequence count for 'dataset_name'
	extimated_data_list - list of estimated abundance (theta) for 'dataset_name'
	"""
	#Development print
	print "\n\t *** Starting Boxplot Function *** \n\n"
	#normalizing data (Z-score)
	normalized_real_data_list = stats.zscore(real_data_list)
	normalized_extimated_data_list = stats.zscore(extimated_data_list)
	#Development print
	print "\t Real Data: ", real_data_list[0:3], " ... ", real_data_list[-1], ". n = ", len(real_data_list), "."
	print "\t Real Data - Normalized (Z-score): ", normalized_real_data_list[0:3], " ... ", normalized_real_data_list[-1], ". n = ", len(normalized_real_data_list), ".\n"
	print "\t Extimated Data: ", extimated_data_list[0:3], " ... ", extimated_data_list[-1], ". n = ", len(extimated_data_list), "."
	print "\t Extimated Data - Normalized (Z-score): ", normalized_extimated_data_list[0:3], " ... ", normalized_extimated_data_list[-1], ". n = ", len(normalized_extimated_data_list), ".\n"
	#Set-up data
	data = [normalized_real_data_list, normalized_extimated_data_list]
	#Creating Boxplot
	plt.figure()
	ax = plt.axes()
	plt.hold(True)
	bp = plt.boxplot(data)
	#Set-up boxplot apparence
	plt.setp(bp['boxes'][0], color='red')
	plt.setp(bp['caps'][0], color='red')
	plt.setp(bp['caps'][1], color='red')
	plt.setp(bp['whiskers'][0], color='red')
	plt.setp(bp['whiskers'][1], color='red')
	plt.setp(bp['fliers'][0], color='red')
	plt.setp(bp['fliers'][1], color='red')
	plt.setp(bp['medians'][0], color='red')
	plt.setp(bp['boxes'][1], color='green')
	plt.setp(bp['caps'][2], color='green')
	plt.setp(bp['caps'][3], color='green')
	plt.setp(bp['whiskers'][2], color='green')
	plt.setp(bp['whiskers'][3], color='green')
	plt.setp(bp['fliers'][2], color='green')
	plt.setp(bp['fliers'][3], color='green')
	plt.setp(bp['medians'][1], color='green')
	ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
	ax.set_axisbelow(True)
	ax.set_title('Boxplot comparison - {0}\n- n={1} redundant IS -'.format(dataset_name, len(real_data_list)))
	#ax.set_xlabel('Datasets')
	ax.set_ylabel('Abundance / SeqCount (Z-score normalized)')
	# draw temporary green and red lines and use them to create a legend
	hR, = plt.plot([1,1],'r-') #extimated
	hG, = plt.plot([1,1],'g-') #real
	plt.legend((hR, hG),('1) Extimated Data - Abundance', '2) Real Data - SequenceCount'))
	hR.set_visible(False)
	hG.set_visible(False)
	#save figure and #show
	plt.savefig(dataset_name + '.comparative_boxplot' + '.pdf', format='pdf')
	#Return
	return 0


def expected_lengths_given_theta (length_phi, freq_phi, theta):
	expected_lenght = 0
	i=0
	for leng in length_phi:
		expected_lenght = expected_lenght + (1 - np.exp(-1*theta*freq_phi[i]))
		i+=1
	return expected_lenght


def phi_VS_theta (length_phi, freq_phi,dataset_name):
	# Retrieving points to plot
	x_expected_unique_lenghts = []
	y_number_of_parent_fragments = []
	N = 1000
	for theta in np.linspace(0, 10000, N):
		y_number_of_parent_fragments.append(theta)
		x_expected_unique_lenghts.append(expected_lengths_given_theta(length_phi, freq_phi, theta))
	#Plot
	plt.figure()
	plt.plot(x_expected_unique_lenghts, y_number_of_parent_fragments, 'b', hold=True)
	plt.xlabel('expected unique lengths (phi-i)')
	plt.ylabel('number of parent fragments (theta-i)')
	plt.title('phi VS theta')
	plt.savefig(dataset_name + '.phiVStheta' + '.pdf', format='pdf')
	return 0




#########################################################################################
### MAIN
#########################################################################################

def main():
	"""
	Main part of the program.
	"""
	#first check args and file paths
	checkArgs(args)
	
	data = args.dataset_file
	f_name = data.split(".")
	print "\n[AP]\t"+"######## "+f_name[0] + '.' + f_name[1]+" ########"
	print "\n[AP]\tChecked inputs, now acquiring data"

	locations_list, length_list = generateDataset(data)

	# Alias for estAbund calling
	estAbund = sonicLength.estAbund

	# Call estAbund and store returned object in results
	results = estAbund(robjects.StrVector(locations_list), robjects.FloatVector(length_list))

	# dev print
	#print (results.r_repr())

	# Put estimation for theta in estimations_theta and associated locations in locations_theta
	theta = results.rx2("theta")
	estimations_theta = tuple(theta)
	locations_theta = tuple(theta.names)

	# Put different fragment lengths in length_phi and associated frequencies in freq_phi
	phi = results.rx2("phi")
	freq_phi = tuple(phi)
	length_phi = tuple(phi.names)

	nameFile = data[0:-20]
	length_phi_numbers = fragmentsLengthPlot(length_phi,freq_phi,length_list,nameFile)

	printThetaInfo(estimations_theta,locations_theta,nameFile)

	host = "172.25.39.2"
	user = "readonly"
	passwd = "readonlypswd"
	db = "sequence_qlam"
	db_table = "osr_p16"
	destfile = nameFile + ".sequence_count" + ".tsv"
	sequence_count = querySeqCount(host,user,passwd,db,db_table,destfile,nameFile)

	sequence_count_list = []
	for v in sequence_count.values():
		sequence_count_list.append(int(v))

	box_plot(sequence_count_list, estimations_theta, nameFile)

	phi_VS_theta(length_phi, freq_phi, nameFile)

	print "\n[AP]\tTask Finished, closing.\n"


# sentinel
if __name__ == "__main__":
    main()