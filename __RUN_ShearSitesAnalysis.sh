#!/bin/bash

echo "
  +--------------------------------------------------------+
  |                                                        |
  |                   Shear Site Analysis                  |
  |                                                        |
  +--------------------------------------------------------+
  |  Author:   Giulio Spinozzi                             |
  |  Date:     October 2013                                |
  |  Version:  1.0                                         |  
  |  Contact:  spinozzi.giulio@hsr.it                      |
  +--------------------------------------------------------+

  REQUIRED VARS and relative ORDER POSITION
        1. ShearSites_identification.py 
            (Python script to generate the Shear Sites for each couple of barcode)   
        2. LTRXX.LCXX.shearsites.tsv 
            (the output file of the previous python script with R1-R2 data of Shear Sites)
        3. abundance_estimation.py 
            (Python script to generate abundance estimation with Berry's model in R)
"


##### ============================ RUN INFO ============================== #####
RUN_STARTED_AT=`date +"%Y-%m-%d %H:%M:%S"`;
RUN_ID="`whoami`"" ${RUN_STARTED_AT}";
#==============================================================================#



##### ============== PYTHON: ShearSites_identification.py ================ #####
echo "PYTHON: ShearSites Identification"
for k in $(ls *.sorted.bed); do 
  n=${k:0:10};
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
  python abundance_estimation.py --dataset $k;
done
#==============================================================================#