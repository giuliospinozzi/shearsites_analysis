#!/bin/bash

echo "
  +----------------------------------------------------------+
  |                                                          |
  |                   Shear Site Analysis                    |
  |                                                          |
  +----------------------------------------------------------+
  |  Author:   Giulio Spinozzi, Brasca Stefano               |
  |  Date:     August 2015                                   |
  |  Version:  1.0                                           |  
  |  Contact:  spinozzi.giulio@hsr.it, brasca.stefano@hsr.it |
  +----------------------------------------------------------+

  REQUIRED FILES:
        - ShearSites_identification.py 
            (Python script to generate the Shear Sites
             for each couple of barcode)   
        - create_matrix command available
            (Integration_Analysis.py)
        - ShearSites_ISaggregation.py
            (Python script to apply fragment mapping aggregation
             around ISs)
        - ShearSites_lengthCorrection.py
            (Python script to apply fragment lengths correction
             after aggregation around ISs)

  REQUIRED VARS and relative ORDER POSITION:
        1. db target schema
        2. db target table
"


##### ============================ RUN INFO ============================== #####
RUN_STARTED_AT=`date +"%Y-%m-%d %H:%M:%S"`;
RUN_ID="`whoami`"" ${RUN_STARTED_AT}";
TODAY=`date +"%Y%m%d%H%M"`;
#==============================================================================#


echo "<< ${RUN_ID} >>
"

##### ========================== INPUT PARAMETERS ======================== #####
DISEASE="${1}";   # AssayValidation
PATIENT="${2}";   # CEMJY
POOL="${3}";      # LANE_1
DBSCHEMA="${4}";  # sequence_assayvalidation
DBTABLE="${5}";   # LANE_1
FASTQDIR="${6}";  # /storage/dx/ngs/data/ShearSites/data/20150706_GSK_HiSeq_validation/ftp.hsr.it/HSR_TIGET/AssayValidation
R2_FASTQ="${7}";  # lane1_NoIndex_L001_R2.fastq.gz
#==============================================================================#


##### ========================== FIXED PARAMETERS ======================== #####
BASEDIR="/opt/NGS/results/$DISEASE/$PATIENT";
OUTDIR="$BASEDIR/quantification/$POOL";
QF="100"
#==============================================================================#


usage="$(basename "$0") <disease> <patient> <pool> <db_schema> <db_table> <fastq_dir> <r2_fastq_file>

where:
    <disease>        Disease ID
    <patient>        Patient ID
    <pool>           Pool
    <db_schema>      DB Schema
    <db_table>       DB Table
    <fastq_dir>      Path to FastQ files
    <r2_fastq_file>  R2 FastQ file (gzipped)
    "

echo "$usage"


##### ========================= PRELIMINARY CHECKS ======================= #####
## Space
LEFTSPACE=`df | grep /dev/sdb1 | awk '{print $4}'`
if [ ${LEFTSPACE} -lt 10000000 ]; then 
    echo "
    
    *****************************************************************
    |                                                               |
    |   YOU DO NOT HAVE ENOUGH SPACE LEFT ON HARD DISK!!            |
    |                                                               |
    |   I WILL NOT PROCEED...                                       |
    |                                                               |       
    |   FREE SPACE BEFORE, LEAVING AT LEAST 10GB ON MAIN DISK       |
    |                                                               |
    *****************************************************************
    
        "; 
    exit;
fi

## Arguments
if [ -z "$1" ]
  then
    echo "----------------> [ERROR] No DISEASE Supplied <----------------
    ";
    exit;
fi
if [ -z "$2" ]
  then
    echo "----------------> [ERROR] No PATIENT Supplied <----------------
    ";
    exit;
fi
if [ -z "$3" ]
  then
    echo "------------------> [ERROR] No POOL Supplied <-----------------
    ";
    exit;
fi
if [ -z "$4" ]
  then
    echo "----------------> [ERROR] No DBSCHEMA Supplied <---------------
    ";
    exit;
fi
if [ -z "$5" ]
  then
    echo "-----------------> [ERROR] No DBTABLE Supplied <---------------
    ";
    exit;
fi
if [ -z "$6" ]
  then
    echo "-----------------> [ERROR] No FASTQDIR Supplied <--------------
    ";
    exit;
fi
if [ -z "$7" ]
  then
    echo "-----------------> [ERROR] No R2_FASTQ Supplied <--------------
    ";
    exit;
fi
#==============================================================================#



echo "

---------------------------------------------------------------------------------
                    STARTING PROCESSING AT: $RUN_STARTED_AT
---------------------------------------------------------------------------------
    " 


# With chrM removal
echo "
+-----------------------------------------------------------------+"
echo "PYTHON: ShearSites Identification"
##### ============== PYTHON: ShearSites_identification.py ================ #####
for k in $(ls ${BASEDIR}/bed/$POOL/*.sorted.md.rel.pg.bed); do
  FILENAME=`basename $k`;
  BARCODE=${FILENAME:0:-21};
  grep -v "chrM" $k > $BARCODE.NOchrM.bed;
  grep -v "chrM" ${k:0:-21}.sorted.allr2reads.bed > $BARCODE.sorted.allr2reads.NOchrM.bed;
  python /opt/applications/scripts/shearsites_analysis/ShearSites_identification.py --bed1 $BARCODE.NOchrM.bed --bed2 $BARCODE.sorted.allr2reads.NOchrM.bed;
done
# out files: *.shearsites.tsv
#==============================================================================#

# # All data
# echo "
# +-----------------------------------------------------------------+"
# echo "PYTHON: ShearSites Identification"
# ##### ============== PYTHON: ShearSites_identification.py ================ #####
# for k in $(ls ${BASEDIR}/bed/$POOL/*.sorted.md.rel.pg.bed); do
#   #FILENAME=`basename $k`;
#   #BARCODE=${FILENAME:0:-21};
#   python /opt/applications/scripts/shearsites_analysis/ShearSites_identification.py --bed1 $k --bed2 ${k:0:-21}.sorted.allr2reads.bed;
# done
# # out files: *.shearsites.tsv
# #==============================================================================#

echo "
+-----------------------------------------------------------------+"
echo "FILTER R2-FASTQ BY RANDOM-BARCODE QUALITY (QF >= ${QF})"
##### ================ FILTER R2-FASTQ BY RANDOM BARCODE QUALITY (QF >= ${QF}) ==================== #####
zcat ${FASTQDIR}/${R2_FASTQ} | trimmomatic SE -phred33 "/dev/stdin" "/dev/stdout" CROP:12 | fastq_quality_filter -q 28 -p ${QF} -Q 33 | fastq_to_fasta -Q 33 -n | fasta2csv | cut -d' ' -f1 | gzip > r2.quality.filtered.list.gz
zcat ${FASTQDIR}/${R2_FASTQ} | fqextract_pureheader <(zcat r2.quality.filtered.list.gz) | gzip -c > r2.qualityFiltered.fastq.gz
#=======================================================================================================#

echo "
+---------------------------------------------------------------------------------------+"
echo "BASH: shrink data according to filtered R2-FASTQ and collect final data headers"
##### ================ BASH: shrink data according to filtered R2-FASTQ and collect final data headers ==================== #####
for k in $(ls *.shearsites.tsv); do
  FILENAME=`basename $k`;
  BARCODE=${FILENAME:0:-15};
  cat $k | extract_header <(zcat r2.quality.filtered.list.gz) > $BARCODE.shearsites.qf${QF}.tsv
  if [ -s $BARCODE.shearsites.qf${QF}.tsv ]
    then 
      awk '{print $1}' $BARCODE.shearsites.qf${QF}.tsv > $BARCODE.headers.txt;
      awk '{print $1}' $BARCODE.shearsites.qf${QF}.tsv >> all_headers.txt;
    else
      rm $BARCODE.shearsites.qf${QF}.tsv;
      echo -e -n "\t > $BARCODE.shearsites.qf${QF}.tsv removed beacuse empty after QF."
  fi
done
#==============================================================================================================================#

echo "
+-----------------------------------------------------------------+"
echo "FQEXTRACT: shrink filtered R2-FASTQ pairing with data"
##### ================ FQEXTRACT: shrink filtered R2-FASTQ pairing with data ==================== #####
zcat r2.qualityFiltered.fastq.gz | fqextract_pureheader all_headers.txt | gzip -c > R2_selected_reads.fastq.gz;
echo "
Removing temp files..."
rm r2.quality.filtered.list.gz;  # rm tmp list of r2 reads to keep by quality
rm r2.qualityFiltered.fastq.gz;  # rm tmp fastq of r2 reads to keep by quality, now replaced by R2_selected_reads.fastq.gz
echo "...done!"
#=====================================================================================================#

echo "
+-----------------------------------------------------------------+"
echo "PYTHON: Launch create_matrix"
##### =================== PYTHON: Launch create_matrix ==================== #####
create_matrix --dbDataset "$DBSCHEMA.$DBTABLE,sequence_qlam.cem_reference" --columns sample,tissue,treatment,vector,enzyme --IS_method gauss --interaction_limit 4 --alpha 0.3 --statistics --tsv --collision;
#===============================================================================#

echo "
+-----------------------------------------------------------------+"
echo "BASH: extract ISs bed data"
##### ================ BASH: extract ISs from matrix ================= #####
awk '{print "chr"$3"\t"$5"\t"$6"\t\t\t"$4}' *${DBTABLE}_StatREPORT_ISs-File*.tsv | tail -n +2 > ISrange.bed;
awk '{print "chr"$3"\t"$9"\t"$9"\t\t\t"$4}' *${DBTABLE}_StatREPORT_ISs-File*.tsv | tail -n +2 > IScoordinate.bed;
#==========================================================================#

echo "
+-----------------------------------------------------------------+"
echo "PYTHON: apply IS aggregation to ShearSite data"
##### ========= PYTHON: apply IS aggregation to ShearSite data ============= #####
# input files: *.shearsites.tsv + ISrange.bed,IScoordinate.bed  (or set args)
python ShearSites_ISaggregation.py --inFiles_nameEnd ".shearsites.qf${QF}.tsv" --outFile_nameEnd ".shearsites.ISfixed.qf${QF}.tsv";  
# out filename: *.shearsites.ISfixed.tsv (or set args)
#=================================================================================#

echo "
+-----------------------------------------------------------------+"
echo "PYTHON: apply length correction to ShearSite data"
##### ========== PYTHON: apply length correction to ShearSite data =========== #####
# input files: *.shearsites.tsv, *.shearsites.ISfixed.tsv (or set args)
python ShearSites_lengthCorrection.py --inFiles_nameEnd ".shearsites.ISfixed.qf${QF}.tsv" --outFile_nameEnd ".shearsites.ISfixed.LENGTHfixed.qf${QF}.tsv";
# out filename: .shearsites.ISfixed.LENGTHfixed.tsv (or set args)
#=================================================================================#

echo "
+-------------------------------------------------------------------------+"
echo "TRIMMOMATIC: extract random barcodes from selected R2 fastq file"
##### ========== TRIMMOMATIC: extract random barcodes from selected R2 fastq file ============ #####
for k in $(ls *.headers.txt); do
  FILENAME=`basename $k`;
  BARCODE=${FILENAME:0:-12};
  echo -e -n "\t > processing ${BARCODE} ... ";
  zcat R2_selected_reads.fastq.gz | fqextract_pureheader $FILENAME | trimmomatic SE -phred33 "/dev/stdin" "/dev/stdout" HEADCROP:12 | fastq_to_fasta -n -Q33 -o ${BARCODE}.randomBC.qf${QF}.fasta;
  echo "done.";
done
echo "
Removing temp files..."
rm R2_selected_reads.fastq.gz;  # not needed anymore
echo "...done!"
#==========================================================================================#

echo "
+-----------------------------------------------------------------+"
echo "PYTHON: build ShearSites data with random barcodes"
##### ========== PYTHON: build shear ShearSites data with random barcodes ============================================================================= #####
# input files: *.shearsites.ISfixed.LENGTHfixed.tsv, *.randomBC.fasta
#              (or set args, also compatible with all *.shearsites.[...].tsv files)
python ShearSites_RandomBC_buildData.py --origFile_nameEnd ".shearsites.ISfixed.LENGTHfixed.qf${QF}.tsv" --BCfasta_nameEnd ".randomBC.qf${QF}.fasta" --outFile_nameEnd ".shearsites.ISfixed.LENGTHfixed.randomBC.qf${QF}.tsv";
# out filename: *.shearsites.ISfixed.LENGTHfixed.randomBC.tsv
#              (or set args properly looking input files names)
#============================================================================================================================================================#

echo "
**********************************************************"
echo "Moving output files in the destination directory..."

mkdir $BASEDIR/quantification
mkdir ${OUTDIR}
mkdir ${OUTDIR}/RandomBC
cp *.pdf *.tsv *.xlsx *.log ${OUTDIR}/RandomBC
echo "...done!"
echo "
Cleaning-up temp files left..."
rm *.pdf *.tsv *.bed *.xlsx *.txt *.randomBC.qf${QF}.fasta
echo "...Finished!
**********************************************************"

echo "

---------------------------------------------------------------------------------
                      ENDING PROCESSING AT: `date +'%Y-%m-%d %H:%M:%S'`
---------------------------------------------------------------------------------
     "
