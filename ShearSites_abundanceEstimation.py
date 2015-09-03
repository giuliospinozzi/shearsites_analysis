#!/usr/bin/python

import argparse, os, sys

import rpy2.robjects as robjects
import numpy as np
import math
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt

# Import sonicLength package (Rpy2 2.1+)
from rpy2.robjects.packages import importr
sonicLength = importr("sonicLength")



header = """

+--------------------------------------------------------------------+
 Author:    Stefano Brasca, Giulio Spinozzi
 Date:    August 2015
 Contact:    brasca.stefano, spinozzi.giulio @hsr.it
 Version:    1.0
+--------------------------------------------------------------------+

 Description
  - This wrapper implements R SonicLength package
  
 Note:
  - The script needs the dataset file
      LTRXX.LCYY.shearsites.(fixed.)sort.uniq.txt,
    and a valid db_schema, db_table

 Steps
    1. Loads the input file
    2. Calls the Berry's Model in R
    3. Saves the output of the model in Python structures
    4. Gives the output files with Phi and Theta
""" 

description = "This is a wrapper with rpy2 for simulate the Berry's model in R."

usage_example = """
Examples of usage:
APP --dataset LTRXX.LCYY.shearsites.(fixed.)sort.uniq.txt --db_schema sequence_something --db_table a_table
"""

# print header

parser = argparse.ArgumentParser(usage = usage_example, epilog = "[ hSR-TIGET - Vector Integration Core - Bioinformatics ] \n", description = description)
parser.add_argument('--dataset', dest="dataset_file", help="Dataset file to process in R. No default option.", action="store", required=True)
parser.add_argument('--db_schema', dest="db_schema", help="Database schema. No default option.", action="store", required=True)
parser.add_argument('--db_table', dest="db_table", help="Database table. No default option.", action="store", required=True)
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


def generateDataset(dataset_file):
    # Load dataset
    file_object = open(dataset_file, 'r')
    dataset = file_object.readlines()
    file_object.close()
    locations_list = []
    length_list = []
    for row in dataset:
        row = row.rstrip('\n')
        row_split = row.split('\t')
        locations_list.append(row_split[0])
        length_list.append(row_split[1])
    return locations_list, length_list


def fragmentsLengthPlot(length_phi,freq_phi,length_list,TAG,dataset_label):
    length_phi_numbers = []
    for num in length_phi:
        length_phi_numbers.append(float(num))
    plt.figure()
    # Phi plot
    plt.plot(length_phi_numbers, freq_phi, 'r', hold=True, label="estimated distribution")
    length_list_numbers = []
    for num in length_list:
        length_list_numbers.append(float(num))
    binning = int(math.ceil(len(set(length_list_numbers))/2))
    # Plot (length-frequency) of input data
    plt.hist(length_list_numbers, bins=binning, normed=True, facecolor='green', alpha=0.3, hold=True, label="real distribution - histogram")
    density = gaussian_kde(length_list_numbers)
    xs = np.linspace(0,max(length_list_numbers)+25,len(set(length_list_numbers))*10+250)
    # Gaussian kde plot
    plt.plot(xs,density(xs), hold=True, label="real distribution - gaussian kde")
    fileName = dataset_label + "." + TAG + ".fragmentsLengthPlot.pdf"
    plt.legend()
    # Labels
    plt.xlabel('fragments length')
    plt.ylabel('probability')
    plt.title(dataset_label + "." + TAG + ' Fragments length data')
    plt.savefig(fileName, format="pdf")
    return length_phi_numbers


def printThetaInfo(estimations_theta,locations_theta,TAG):
    f_out = open(TAG+".theta.tsv","w")
    dict_theta = {}
    for l in range(0,len(locations_theta)):
        dict_theta[locations_theta[l]] = estimations_theta[l]
        f_out.write(str(locations_theta[l]) + '\t' + str(estimations_theta[l]) + '\n')


def queryDataset(host,user,passwd,db,db_table,destfile,TAG):
    # System Call
    query = "mysql -h %(host)s -u %(user)s --password=%(passwd)s %(db)s --skip-column-names -e \"SELECT DISTINCT sample, tissue, treatment, vector, enzyme FROM %(db_table)s WHERE tag='%(TAG)s'\" > %(destfile)s" %{
     'host': host,
     'user': user,
     'passwd': passwd,
     'db': db,
     'db_table': db_table,
     'TAG': TAG,
     'destfile': destfile,
    }
    os.system(query)
    f_in = open(destfile, 'r')
    value = None
    for line in f_in:
        string_splitted = line.split('\t')
        value = str(string_splitted[0]) + "." + str(string_splitted[1]) + "." + str(string_splitted[2]) + "." + str(string_splitted[3]) + "." + str(string_splitted[4])
    f_in.close()
    os.remove(destfile)
    return value


def expected_lengths_given_theta (length_phi, freq_phi, theta):
    expected_lenght = 0
    i=0
    for leng in length_phi:
        expected_lenght = expected_lenght + (1 - np.exp(-1*theta*freq_phi[i]))
        i+=1
    return expected_lenght


def phi_VS_theta (length_phi, freq_phi,TAG,dataset_label):
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
    plt.savefig(dataset_label + "." + TAG + '.phiVStheta' + '.pdf', format='pdf')
    return 0


#def is_CEM (genome_location_string):
#    response = False
#    cem_symbol = None
#    cem_locations = [('chr11',64536791,64537194), ('chr17',47732128,47732367), ('chr17',2032108,2032381), ('chr2',24546189,24546598), ('chr2',73761972,73762443), ('chr16',28497050,28497538)]
#    #cem_locations.append(('chr1',1776600,1776700))    #Fake to test: 1776639
#    cem_symbols_list = ['$','#','@','*','%','!']
#    
#    splitted_location = genome_location_string.split(' ') # 0 -> chr; 1 -> locus; 2 -> strand
#    splitted_location[1] = int(splitted_location[1])
#
#    i=0;
#    for cem_location_tupla in cem_locations:
#        if (splitted_location[0] in cem_location_tupla[0]):
#            if (splitted_location[1] in range(cem_location_tupla[1],cem_location_tupla[2]+1)):
#                response = True
#                cem_symbol = cem_symbols_list[i]
#                return response, cem_symbol
#        i+=1
#
#    return response, cem_symbol




#########################################################################################
### MAIN
#########################################################################################

def main():
    """
    Main part of the program.
    """
    # Check args and file paths
    checkArgs(args)
    
    # Screen print for user
    dataset_file = args.dataset_file
    dataset_file_split = dataset_file.split(".")
    print "\n[AP]\t" + "######## " + dataset_file_split[0] + '.' + dataset_file_split[1] + " ########"
    print "\n[AP]\tInput checked, now acquiring data ... "
    
    # DB connection parameter
    host = "localhost"
    user = "readonly"
    passwd = "readonlypswd"
    db = args.db_schema
    db_table = args.db_table
    
    # Get human readable dataset_label from TAG, through DB
    TAG = dataset_file[0:-31]  #remove '.shearsites.fixed.sort.uniq.txt'
    dataset_label = queryDataset(host, user, passwd, db, db_table, "{TAG}_dataset_label.DB.tmp".format(TAG=str(TAG)), TAG)
    
    if dataset_label is not None:
        dataset_label = dataset_label.rstrip('\n')
        dataset_label = dataset_label.replace("/","-")
        
        locations_list, length_list = generateDataset(dataset_file)
        
        if len(locations_list) < 2:
            print "\n[AP]\t{dataset_label} has only one unique line! Can't estimate anything.".format(dataset_label=str(dataset_label))
            print "[AP]\tSKIP THIS DATASET!\n"
            return 0

        # Alias for estAbund calling
        estAbund = sonicLength.estAbund

        # Call estAbund and store returned object in results
        results = estAbund(robjects.StrVector(locations_list), robjects.FloatVector(length_list))

        # Put estimation for theta in estimations_theta and associated locations in locations_theta; then organize data in dic_of_theta
        theta = results.rx2("theta")
        estimations_theta = tuple(theta)
        locations_theta = tuple(theta.names)
        # dic_of_theta
        dic_of_theta = {}
        for i in range(len(locations_theta)):
            dic_of_theta.update({locations_theta[i]:estimations_theta[i]})

        # Put different fragment lengths in length_phi and associated frequencies in freq_phi
        phi = results.rx2("phi")
        freq_phi = tuple(phi)
        length_phi = tuple(phi.names)

        length_phi_numbers = fragmentsLengthPlot(length_phi,freq_phi,length_list,TAG,dataset_label)

        printThetaInfo(estimations_theta,locations_theta,TAG)

        # Plot: unique lengths retrieved for a genomic location VS expected number of parent fragment for the same location
        phi_VS_theta(length_phi, freq_phi, TAG, dataset_label)

        # Last screen print for user
        print "\n[AP]\tTask Finished, closing.\n"
    else:
        # This is usually an error, print some debug info
        print "\n[AP]\tCan't find metadata in DB!"
        print "\t\tDataset details:"
        print "\t\t* dataset_file: {dataset_file}".format(dataset_file=str(dataset_file))
        print "\t\t* extracted TAG: {TAG}".format(TAG=str(TAG))
        print "\t\tDB details:"
        print "\t\t* host: {host}".format(host=str(host))
        print "\t\t* user: {user}".format(user=str(user))
        print "\t\t* passwd: {passwd}".format(passwd=str(passwd))
        print "\t\t* db: {db}".format(db=str(db))
        print "\t\t* db_table: {db_table}".format(db_table=str(db_table))
        print "\t\tDataset Label returned by query:"
        print "\t\t* dataset_label: {dataset_label}".format(dataset_label=str(dataset_label))
        print "\n[AP]\tSKIP THIS DATASET!\n"
        
    return 0


# sentinel
if __name__ == "__main__":
    main()