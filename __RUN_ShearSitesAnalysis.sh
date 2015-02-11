#!/bin/bash

echo "
  +--------------------------------------------------------+
  |                                                        |
  |                   Shear Site Analysis                  |
  |                                                        |
  +--------------------------------------------------------+
  |  Author:   Giulio Spinozzi                             |
  |  Date:     January 2015                                |
  |  Version:  3.0                                         |  
  |  Contact:  spinozzi.giulio@hsr.it                      |
  +--------------------------------------------------------+

  REQUIRED FILES:
        - ShearSites_identification.py 
            (Python script to generate the Shear Sites for each couple of barcode)   
        - LTRXX.LCXX.shearsites.tsv 
            (the output file of the previous python script with R1-R2 data of Shear Sites)
        - abundance_estimation.py 
            (Python script to generate abundance estimation with Berry's model in R)

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
OUTDIR="$BASEDIR/quantification";
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


##### ============== PYTHON: ShearSites_identification.py ================ #####
echo "
+--------------------------------------------------------+"
echo "PYTHON: ShearSites Identification"
for k in $(ls ${BASEDIR}/bed/$POOL/*.sorted.md.rel.pg.bed); do
  n=${k:0:-21};
  python /opt/applications/scripts/shearsites_analysis/ShearSites_identification.py --bed1 $k --bed2 $n.sorted.allr2reads.bed;
done
#==============================================================================#


echo "
+--------------------------------------------------------+"
echo "BASH: output txt files for R model"
##### =============== BASH: output txt files for R model ================= #####
for k in $(ls *.shearsites.tsv); do 
  awk '{print $2" "$3" "$6"\t"$12}' $k | sort | uniq > ${k:0:-4}.uniq.txt;
done

for k in $(ls *.shearsites.tsv); do 
  awk '{print $2" "$3" "$6"\t"$12}' $k | sort > ${k:0:-4}.txt;
done
#==============================================================================#


echo "
+--------------------------------------------------------+"
echo "RPY2: Abundance Estimation with Berry's Model in R"
##### ================ PYTHON: abundance_estimation.py =================== #####
for k in $(ls *.shearsites.uniq.txt); do 
  python /opt/applications/scripts/shearsites_analysis/abundance_estimation.py --dataset $k --db_schema $DBSCHEMA --db_table $DBTABLE;
done
#==============================================================================#


echo "
**********************************************************"
echo "Moving output files in the destination directory..."
mkdir ${OUTDIR}
cp *.tsv *.txt *.pdf ${OUTDIR}
echo "...done!"
echo "
Removing temp files..."
rm *.pdf *.txt *.tsv
echo "...Finished!
**********************************************************"


echo "

---------------------------------------------------------------------------------
                      ENDING PROCESSING AT: `date +'%Y-%m-%d %H:%M:%S'`
---------------------------------------------------------------------------------
    " 