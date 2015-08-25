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
        - LTRXX.LCXX.shearsites.tsv 
            (the output file of the previous python script
             with R1-R2 data of Shear Sites)
		- create_matrix command available
			(Integration_Analysis.py)
		- ShearSites_lengthCorrection.py
			(Python script to apply fragment lengths correction
			  after grouping around ISs)
        - ShearSites_abundanceEstimation.py
            (Python script to estimate abundance through
              sonicLength R package)

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
DISEASE="${1}";  #MLD
PATIENT="${2}";  #MLD01
POOL="${3}";     #pool1
DBSCHEMA="${4}"; #sequence_qlam
DBTABLE="${5}";  #new_method
#==============================================================================#


##### ========================== FIXED PARAMETERS ======================== #####
BASEDIR="/opt/NGS/results/$DISEASE/$PATIENT";
OUTDIR="$BASEDIR/quantification/$POOL";
BEDTOOLS="/opt/applications/bin/bedtools/bedtools-2.17.0/bin/bedtools";
#==============================================================================#


usage="$(basename "$0") <disease> <patient> <pool> <db_schema> <db_table>

where:
    <disease>    Disease ID
    <patient>    Patient ID
    <pool>       Pool
    <db_schema>  DB Schema
    <db_table>   DB Table
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
==============================================================================#


echo "

---------------------------------------------------------------------------------
                    STARTING PROCESSING AT: $RUN_STARTED_AT
---------------------------------------------------------------------------------
    " 


echo "
+--------------------------------------------------------+"
echo "PYTHON: ShearSites Identification"
##### ============== PYTHON: ShearSites_identification.py ================ #####
for k in $(ls ${BASEDIR}/bed/$POOL/*.sorted.md.rel.pg.bed); do
  n=${k:0:-21};
  python /opt/applications/scripts/shearsites_analysis/ShearSites_identification.py --bed1 $k --bed2 $n.sorted.allr2reads.bed;
done
# out files: *.shearsites.tsv
#==============================================================================#

echo "
+--------------------------------------------------------+"
echo "BASH: extract ShearSite bed data"
##### ================ BASH: extract ShearSite bed data =================== #####
for k in $(ls *.shearsites.tsv); do
  awk '{print $2"\t"$3"\t"$3"\t\t\t"$6"\t"$12}' $k | sort | sort -t $'\t' -k1,1 -k2,2 | uniq > ${k:0:-4}.raw.sort.uniq.bed;
done
# out files: *.shearsites.raw.sort.uniq.bed
#==============================================================================#



echo "
+--------------------------------------------------------+"
echo "PYTHON: Launch create_matrix"
##### =================== PYTHON: Launch create_matrix ==================== #####
create_matrix --dbDataset "$DBSCHEMA.$DBTABLE,sequence_qlam.cem_reference" --columns sample,tissue,treatment,vector,enzyme --IS_method gauss --interaction_limit 4 --alpha 0.3 --statistics --tsv --collision
#===============================================================================#

echo "
+--------------------------------------------------------+"
echo "BASH: extract ISs coordinates bed data"
##### ========= BASH: extract ISs coordinates from matrix =========== #####
awk '{print "chr"$3"\t"$5"\t"$6"\t\t\t"$4}' *${DBTABLE}_StatREPORT_ISs-File*.tsv | tail -n +2 > ISrange.bed
awk '{print "chr"$3"\t"$9"\t"$9"\t\t\t"$4}' *${DBTABLE}_StatREPORT_ISs-File*.tsv | tail -n +2 > IScoordinate.bed
#==========================================================================#



echo "
+--------------------------------------------------------+"
echo "BASH: apply IS correction to ShearSite data"
##### ========= BASH: apply IS correction to ShearSite data ============= #####
for k in $(ls *.shearsites.raw.sort.uniq.bed); do
  # Old version
  # $BEDTOOLS intersect -a $k -b ISrange.bed -wb -s > ${k:0:-4}.tmp.txt;
  # awk '{print $6"\t"$7"\t"$8"\t\t\t"$9"\t"$5}' ${k:0:-4}.tmp.txt > ${k:0:-4}.tmp.bed;
  # $BEDTOOLS intersect -a ${k:0:-4}.tmp.bed -b IScoordinate.bed -wb -s > ${k:0:-4}.txt;
  # awk '{print $6"\t"$7"\t"$8"\t\t\t"$9"\t"$5}' ${k:0:-4}.txt | sort | uniq > ${k:0:-29}.shearsites.sort.uniq.bed;
  $BEDTOOLS intersect -a $k -b ISrange.bed -wb -s | awk '{print $6"\t"$7"\t"$8"\t\t\t"$9"\t"$5}' > ${k:0:-4}.tmp.bed;
  $BEDTOOLS intersect -a ${k:0:-4}.tmp.bed -b IScoordinate.bed -wb -s | awk '{print $6"\t"$7"\t"$8"\t\t\t"$9"\t"$5}' > ${k:0:-29}.shearsites.ISfixed.bed;
done
# *.raw.sort.uniq.bed -> *.shearsites.ISfixed.bed
# ...not properly a bed, lengths are put after strand in 7th field (thickStart)
#=============================================================================#

echo "
+--------------------------------------------------------+"
echo "PYTHON: apply length correction to ShearSite data"
##### ========== PYTHON: apply length correction to ShearSite data =========== #####
# input files: *.shearsites.raw.sort.uniq.bed, *.shearsites.ISfixed.bed (or set args)
python ShearSites_lengthCorrection.py
# out filename: *.shearsites.ISfixed.LENGTHfixed.bed (or set args)
#=================================================================================#


echo "
+--------------------------------------------------------+"
echo "BASH: extract data for sonicLength"
##### =============== BASH: extract data for sonicLength ================= #####
for k in $(ls *.shearsites.ISfixed.LENGTHfixed.bed); do
  awk '{print $1"\t"$2"\t"$4"\t"$5}' $k | sort -u -t $'\t' -k1,1 -k2,2 -k3,3 -k4,4n | awk '{print $1" "$2" "$3"\t"$4}' > ${k:0:-35}.shearsites.fixed.sort.uniq.txt;
done
#out filename: *.shearsites.fixed.sort.uniq.txt - sonicLength input
#==============================================================================#


echo "
+--------------------------------------------------------+"
echo "RPY2: Abundance Estimation with Berry's Model in R"
##### ================ PYTHON: abundance_estimation.py =================== #####
for k in $(ls *.shearsites.fixed.sort.uniq.txt); do 
  python ShearSites_abundanceEstimation.py --dataset $k --db_schema $DBSCHEMA --db_table $DBTABLE;
done
#==============================================================================#



echo "
**********************************************************"
echo "Moving output files in the destination directory..."
mkdir $BASEDIR/quantification
mkdir ${OUTDIR}
cp *.tsv *.txt *.pdf ${OUTDIR}
# echo "...done!"
# echo "
# Removing temp files..."
# rm *.pdf *.txt *.tsv *.bed *.xlsx
echo "...Finished!
**********************************************************"



echo "

---------------------------------------------------------------------------------
                      ENDING PROCESSING AT: `date +'%Y-%m-%d %H:%M:%S'`
---------------------------------------------------------------------------------
    " 