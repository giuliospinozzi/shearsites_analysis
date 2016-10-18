# -*- coding: utf-8 -*-


"""
Created on Thu Oct 13 14:48:59 2016

@author: stefano
"""


import sys
import os
import argparse


description = """
+---+ GENERAL DESCRIPTION +---------------------------------------------------+

Given one or more dataset to analyze, this program computes Integration Sites
quantification matrixes.

According to the features of the whole input data, this program yields as
output up to five files:
1) seqCount_matrix.tsv - count of sequencing reads for each IS
2) ShsCount_matrix.tsv - count of distinct shearsites for each IS
3) barcodeCount_matrix.tsv - count of distinct random barcodes for each IS
4) cellCount_matrix.tsv - count of parental fragments of each IS
   (obtained as the count of distinct shearsite&random-barcode couples)
5) fragmentEstimate_matrix.tsv - estimate of parental fragments of each IS
   (computed by sonicLength R-package exploiting shearsite data only)
By default such files are created in the current working directory.

Integration Sites aggregation is performed through a sliding window approach,
centered on the sequence count local maximum. In case of ties, the first locus
hosting a maximum is taken as IS locus (strand-wise).

If available, random barcode's data are filtered by edit-distance trying to
mitigate overestimation biases due to PCR-artifacts / sequencing-erros.

Example of standard usage:
a) {program} /opt/NGS/results,MLD,MLDHE01
   (process all pools found in /opt/NGS/results related to MLD disease and 
   MDLHE01 patient)
b) {program} /opt/NGS/results,Pirazzoli,PE_LungCancerResistance,VP1,VP2
   (process only pools VP1 and VP2 of Pirazzoli PE_LungCancerResistance)
   
Here strings like
a) data_base_dir,disease,patient
b) data_base_dir,disease,patient,a_pool,another_pool,one_pool_more,...
are called 'dataset tuple' and {program} can accept one or more,
separated eachother by a space, to be processed together.

+---+ FURTHER DETAILS +-------------------------------------------------------+

* About Integration Sites aggregation:
Integration Site aggregation of covered bases is performed by default by means
of a 'sliding-window' approach. The maximum gap allowed between two covered
bases in the same window is tunable through '--ISs_max_gap' argument (default:
7) as well as the maximum overall window extent is tunable through
'--ISs_max_span' argument (default: 8). Such rules holds together and thus the
window behaviour is fully determined.
All the covered bases in the same window are merged as they were one, now
called IS: total sequence count is updated, shearsites are corrected, related
distinct random barcodes are recomputed and genomic coordinates are inherited
from the covered base with the highest sequence count; in case of ties, the
first locus hosting a maximum is taken as IS locus (strand-wise).
If you want a simple 'covered base matrix', just use '-cb' / '--covered_bases'
option (note: any given 'IS argument' will be silently ignored).

* About Random Barcode Edit-distance Filter:
Each IS (or covered base, in case of '-cb' / '--covered_bases' option)
undergoes an inspection of the mutual edit distance among its random barcodes. 
This inspection can be done inside its shearsite compartments (default) or
overall (--filter_ignoring_shearsites): the latter may be conceptually wrong
and for sure extremely time&memory-consuming.
Finally, for each comparison yileding and edit-distance equal or greater than 3
(defualt, can be tuned through '--filter_edit_distance_threshold' argument),
the barcode with the lower sequence count is discarded (along with related
data), while the other one (higher sequence count) is kept.
This is done to try to mitigate overestimation biases due to PCR-artifacts /
sequencing-erros, discarding adulterated data.

* About Output:
Output files are plain text file, utf-8, tab-separated. Each row represents an
IS (or covered base, in case of '-cb' / '--covered_bases' option) demultiplexed
in columns according to 'sample IDs' found in data. Such files are created in
the current working directory as default, or in any location passed through
-o /--out_dir_path. The argument must be an absoute path of a directory: if not
exists will be created. Further, -s / --subfolder argument allow to set an
additional subfolder (None by default).

* Other Options:
-q / --quiet option prevent the program from printing pretty anything, except
for main warnings, critical errors and python's exceptions.

* Notes:
Many controls are spreaded inside the code, about:
- python environment
- R environment
- argument parsed settings
- output paths and permissions
- input data compliance
- analysis consistency
However they not aim at being exhaustive. To trust results, be sure about what
are you doing!

+---+ ABOUT ARGUMENTS: +------------------------------------------------------+""".format(program=str(os.path.basename(sys.argv[0])))



class MyParser(argparse.ArgumentParser):
    def error(self, message):
        print ""
        self.print_help()
        print "\n+-----------------------------------------------------------------------------+"
        print "\n[ERROR]: %s" % message
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
        raise argparse.ArgumentTypeError("invalid syntax: {s}. See help.".format(s=str(s)))

def min_limited_int(n, m):
    try:
        n = int(n)
    except:
        raise argparse.ArgumentTypeError("N={n} must be an int.".format(s=str(n)))
    if n < m:
        raise argparse.ArgumentTypeError("N={m} > {n} is the minimum allowed.".format(m=str(m), n=str(n)))
    return int(n)

def max_limited_int(n, M):
    try:
        n = int(n)
    except:
        raise argparse.ArgumentTypeError("N={n} must be an int.".format(s=str(n)))
    if n > M:
        raise argparse.ArgumentTypeError("N={M} < {n} is the maximum allowed'".format(M=str(M), n=str(n)))
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
parser = MyParser(description=description, formatter_class=argparse.RawDescriptionHelpFormatter)
# Set dataset_tuple_list - REQUIRED
parser.add_argument("dataset_tuple_list", metavar='DATA_TUPLES', nargs='+', type=dataset_tuple, help="one (or more) dataset to process. Each dataset is expected to be a tuple composed by at least 3 items, comma separated with no spaces (e.g: data_base_dir,disease,patient); furter item(s) following these three will be interpreted as the only pool(s) to process (e.g.: data_base_dir,disease,patient,poolA,poolB), otherwise all the pools of the 3-items-tuple dataset will be taken. A space separates distinct tuples.")  # return a list
# Set verbose - optional
parser.add_argument("-q", "--quiet", action="store_false", default=True, help="silent execution. (default: verbose).")
# Set common_output_ground_dir - optional
parser.add_argument("-o", "--out_dir_path", metavar='ABS_OUT_PATH', default=os.getcwd(), help="set an absolute directory path where write output matrix files. (default: current working directory).")
# Set matrix_outfolder - optional
parser.add_argument("-s", "--subfolder", metavar='FOLDER_NAME', default='', help="specify a subfolder in ABS_OUT_PATH where write output matrix files. (default: no subfolder).")
### Set do_ISs - optional
parser.add_argument("-cb", "--covered_bases", action="store_false", default=True, help="prevent from ISs aggregation, yielding a simple 'covered bases matrix'. If this argument is given, other '--ISs_...' arguments/settings will be ignored. (default: DO ISs aggregation).")
# Set do_ISs param: ensembles_per_sample - optional
parser.add_argument("--ISs_per_sample", action="store_true", default=False, help="aggregate covered bases in ISs sample-by-sample. (default: aggregate covered bases in ISs across all data).")
# Set ensembles_max_dist - optional
parser.add_argument("--ISs_max_gap", metavar='N', type=positive_int, default=7, help="if (locus[n+1] - locus[n] > N), the IS window is truncated at locus[n]. (default: 7).")
# Set ensembles_max_span - optional
parser.add_argument("--ISs_max_span", metavar='N', type=greater_than_1_int, default=8, help="if (locus[n+1] - locus[0] > N), the IS window is truncated at locus[n]. (default: 8).")
### Set filter_data - optional
parser.add_argument("-uf", "--unfiltered", action="store_false", default=True, help="prevent data from being filtered by random barcode's edit-distance. If this argument is given, other '--filter_...' arguments/settings will be ignored. (default: DO filter by random barcode's edit-distance).")
# Set ED_treshold - optional
parser.add_argument("--filter_edit_distance_threshold", metavar='N', type=range_1_11_int, default=3, help="filter out data whose random barcodes are distant N or less (Levenshtein distance): for each comparison, only the counterpart with the lower sequence count is discarded, while the other one (higher sequence count) is kept.  (default: 3).")
# Set inside_ShS - optional
parser.add_argument("--filter_ignoring_shearsites", action="store_false", default=True, help="apply the filter by edit-distance across all the random barcodes of each IS, even if belonging to different shearsites. (default: apply the filter within shearsite compartments of each IS).")

#Parse
args = parser.parse_args()
# print vars(args)

