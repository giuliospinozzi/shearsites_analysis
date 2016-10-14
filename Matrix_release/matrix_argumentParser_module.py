# -*- coding: utf-8 -*-


"""
Created on Thu Oct 13 14:48:59 2016

@author: stefano
"""


description = """
DESCRIPTION HERE
"""

import sys
import os
import argparse


class MyParser(argparse.ArgumentParser): 
   def error(self, message):
      print 'error: %s\n' % message
      print ""
      self.print_help()
      sys.exit("\n[QUIT]\n")


def dataset_tuple(s, arg_sep=','):
    try:
        arg_split = s.split(arg_sep)
        ground_dir, DISEASE, PATIENT, POOLS = arg_split[0], arg_split[1], arg_split[2], arg_sep.join(arg_split[3:])
        if POOLS == '':
            POOLS = False  # Take all pools if nothing given through 's'
        else:
            POOLS = arg_split[3:]  # List of pools (1 item at least)
        dataset_tuple = ground_dir, DISEASE, PATIENT, POOLS
        return dataset_tuple
    except:
        raise argparse.ArgumentTypeError("ERROR HERE, somthing like 'invalid syntax for {s}'".format(s=str(s)))

def min_limited_int(n, m):
    try:
        n = int(n)
    except:
        raise argparse.ArgumentTypeError("ERROR HERE, somthing like '{n} must be an int'".format(s=str(n)))
    if n < m:
        raise argparse.ArgumentTypeError("ERROR HERE, somthing like '{m} > {n} is the min allowed'".format(m=str(m), n=str(n)))
    return int(n)

def max_limited_int(n, M):
    try:
        n = int(n)
    except:
        raise argparse.ArgumentTypeError("ERROR HERE, somthing like '{n} must be an int'".format(s=str(n)))
    if n > M:
        raise argparse.ArgumentTypeError("ERROR HERE, somthing like '{M} < {n} is the max allowed'".format(M=str(M), n=str(n)))
    return int(n)

def range_limited_int(n, m, M):
    n = min_limited_int(n, m)
    n = max_limited_int(n, M)
    return int(n)

def positive_int(n):
    return min_limited_int(n, 1)
    
def greater_than_1_int(n):
    return min_limited_int(n, 2)
    
def range_1_11_int(n):
    return range_limited_int(n, 1, 11)


#parser = MyParser(description=description, formatter_class=argparse.RawTextHelpFormatter)
parser = MyParser(description=description)
# Set dataset_tuple_list - REQUIRED
parser.add_argument("dataset_tuple_list", metavar='DATA_TUPLES', nargs='+', type=dataset_tuple, help="HERE HELP FOR 'DATA'")  # return a list
# Set verbose - optional
parser.add_argument("-q", "--quiet", action="store_false", default=True, help="silent execution")
# Set common_output_ground_dir - optional
parser.add_argument("-o", "--out_dir_path", metavar='ABS_OUT_PATH', default=os.getcwd(), help="HERE HELP FOR ABS_OUT_PATH")
# Set matrix_outfolder - optional
parser.add_argument("-s", "--subfolder", metavar='FOLDER_NAME', default='', help="HERE HELP FOR FOLDER_NAME")
### Set do_ISs - optional
parser.add_argument("-cb", "--covered_bases", action="store_false", default=True, help="prevent from ISs aggregation, yielding a simple 'covered bases matrix'. In this case, '--ISs_...' arguments will be ignored if given.")
# Set do_ISs param: ensembles_per_sample - optional
parser.add_argument("--ISs_per_sample", action="store_true", default=False, help="aggregate covered bases in ISs sample-by-sample.")
# Set ensembles_max_dist - optional
parser.add_argument("--ISs_max_gap", metavar='N', type=positive_int, default=7, help="if (locus[n+1] - locus[n] > N), the IS window is truncated at locus[n].")
# Set ensembles_max_span - optional
parser.add_argument("--ISs_max_span", metavar='N', type=greater_than_1_int, default=8, help="if (locus[n+1] - locus[0] > N), the IS window is truncated at locus[n].")
### Set filter_data - optional
parser.add_argument("-uf", "--unfiltered", action="store_false", default=True, help="prevent data from being filtered by random barcode's edit-distance.")
# Set ED_treshold - optional
parser.add_argument("--filter_edit_distance_threshold", metavar='N', type=range_1_11_int, default=3, help="apply filter on data whose random barcodes are distant N or less (Levenshtein distance): for each comparison, the counterpart with the higher sequence count is kept. In this case, '--filter_...' arguments will be ignored if given.")
# Set inside_ShS - optional
parser.add_argument("--filter_ignoring_shearsites", action="store_false", default=True, help="prevent the filter by random barcode's edit-distance from gruping integration data by shearsites")

#Parse
args = parser.parse_args()

# print vars(args)












