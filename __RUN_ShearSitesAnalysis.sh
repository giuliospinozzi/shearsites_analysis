#!/bin/bash

echo "
  +--------------------------------------------------------+
  |                                                        |
  |                   Shear Site Analysis                  |
  |                                                        |
  +--------------------------------------------------------+
  |  Author:   Giulio Spinozzi                             |
  |  Date:     October 2013                                |
  |  Version:  2.0                                         |  
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
#==============================================================================#


echo "<< ${RUN_ID} >>
"

##### ============================ SETTINGS ============================== #####
DBSCHEMA="${1}";
DBTABLE="${2}";
PROJECT="ShearSites";
EXPERIMENT="qLam";
NGSWORKINGPATH="/opt/NGS/results";
INPUTDIR_POOL_BED="${NGSWORKINGPATH}/${PROJECT}/${EXPERIMENT}/NASarchive/bed/"${DBTABLE};
OUTDIR="${NGSWORKINGPATH}/${PROJECT}/${EXPERIMENT}/NASarchive/quantification/"${DBTABLE};
#==============================================================================#


##### ============== PYTHON: ShearSites_identification.py ================ #####
echo "
+--------------------------------------------------------+"
echo "PYTHON: ShearSites Identification"
for k in $(ls ${INPUTDIR_POOL_BED}/*.sorted.bed); do
  n=${k:0:-28};
  python ShearSites_identification.py --bed1 $k --bed2 $n.sorted.allr2reads.bed;
done
#==============================================================================#


echo "
+--------------------------------------------------------+"
echo "BASH: output txt files for R model"
##### =============== BASH: output txt files for R model ================= #####
for k in $(ls *.shearsites.tsv); do 
  awk '{print $2" "$3" "$6"\t"$12}' $k | sort | uniq > ${k:0:21}.uniq.txt;
done

for k in $(ls *.shearsites.tsv); do 
  awk '{print $2" "$3" "$6"\t"$12}' $k | sort > ${k:0:21}.txt;
done
#==============================================================================#


echo "
+--------------------------------------------------------+"
echo "RPY2: Abundance Estimation with Berry's Model in R"
##### ================ PYTHON: abundance_estimation.py =================== #####
for k in $(ls *.shearsites.uniq.txt); do 
  n=${k:0:10};
  python abundance_estimation.py --dataset $k --db_schema $DBSCHEMA --db_table $DBTABLE;
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