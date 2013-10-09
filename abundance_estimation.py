#!/usr/bin/python

import argparse, os

# Required
import rpy2.robjects as robjects
import numpy as np
import matplotlib.pyplot as plt

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


def plotPhi(length_phi,freq_phi,nameFig):
	length_phi_numbers = []
	for num in length_phi:
		length_phi_numbers.append(int(num))

	# generate plot
	plt.bar(length_phi_numbers, freq_phi, width=1.0, bottom=None, hold=False)
	fileName = nameFig + "phi.pdf"
   	plt.savefig(fileName, format="pdf")
   	#plt.show()


def plotHistFreq(length_list,nameFig):
	# Plot (length-frequency) of input data
	length_list_numbers = []
	for num in length_list:
		length_list_numbers.append(int(num))
	plt.hist(length_list_numbers, bins=max(length_list_numbers), normed=True, hold=False)
	fileName = nameFig + "histFreq.pdf"
   	plt.savefig(fileName, format="pdf")
   	#plt.show()
	return length_list_numbers


def plotGaussianDensity(length_list_numbers,nameFig):
	# Gaussian Density Plot (length-frequency) of input data
	density = gaussian_kde(length_list_numbers)
	xs = np.linspace(0,400,400)
	plt.plot(xs,density(xs), hold=False)
	fileName = nameFig + "gaussDensity.pdf"
   	plt.savefig(fileName, format="pdf")
   	#plt.show()



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

	data = args.dataset_file

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

	nameFig = data[0:-19]
	plotPhi(length_phi,freq_phi,nameFig)
	length_list_numbers = plotHistFreq(length_list,nameFig)
	plotGaussianDensity(length_list_numbers,nameFig)

	print "\n[AP]\tTask Finished, closing.\n"


# sentinel
if __name__ == "__main__":
    main()