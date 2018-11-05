#!/bin/bash

## ------------------------------------------------------------
## run IEDB Immunogenicity prediction on peptides
## ------------------------------------------------------------

module load biokit

cd ./data/vcf/proteome 

for FILE in `ls -1 *differential*9.fasta`; do
  SAMPLE=`echo $FILE | cut -d "_" -f 1` 
  sed '/^>/d' $FILE > tmp.file
  /proj/2000066/tools/IEDB/predict_immunogenicity.py tmp.file > ./immunogenicity/$SAMPLE.9.immunogenicity
  rm -f tmp.file
done



