#!/bin/bash

# =====================================================
# sbatch array file:
# Run peptide binding affinity for HLA-A & HLA-B
#   
# =====================================================

#SBATCH -J peptide
#SBATCH -o peptide_out_%A_%a.txt
#SBATCH -e peptide_err_%A_%a.txt
#SBATCH -t 48:00:00
#SBATCH --mem-per-cpu=5000
#SBATCH --array=1-161
#SBATCH -n 1
#SBATCH -p serial

module load biokit

cd ./data/vcf/proteome/immunogenicity

FILE=`ls -1 *Imm_Exp_GvHDBM.fasta | sed "${SLURM_ARRAY_TASK_ID}q;d"` 
SAMPLE=`echo $FILE | cut -d '_' -f 1`
HLA_A1=`grep -n $SAMPLE ./data/HLA_types/donor.types | cut -f 2 | cut -d ',' -f 1 | tr -d '[:space:]'`
HLA_A2=`grep -n $SAMPLE ./data/HLA_types/donor.types | cut -f 2 | cut -d ',' -f 2 | tr -d '[:space:]'`
HLA_B1=`grep -n $SAMPLE ./data/HLA_types/donor.types | cut -f 3 | cut -d ',' -f 1 | tr -d '[:space:]'`
HLA_B2=`grep -n $SAMPLE ./data/HLA_types/donor.types | cut -f 3 | cut -d ',' -f 2 | tr -d '[:space:]'`
./src/IEDB/mhc_i/src/predict_binding.py IEDB_recommended $HLA_A1 9 $FILE > $SAMPLE.A1.IEDB.out
./src/IEDB/mhc_i/src/predict_binding.py IEDB_recommended $HLA_A2 9 $FILE > $SAMPLE.A2.IEDB.out
./src/IEDB/mhc_i/src/predict_binding.py IEDB_recommended $HLA_B1 9 $FILE > $SAMPLE.B1.IEDB.out
./src/IEDB/mhc_i/src/predict_binding.py IEDB_recommended $HLA_B2 9 $FILE > $SAMPLE.B2.IEDB.out


