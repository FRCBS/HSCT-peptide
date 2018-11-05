#!/bin/bash

## -----------------------------------------------------------
## 
## Run the scripts for:
## 1. Recipient-mismatched peptide extraction
## 2. Peptide immunogenicity analysis
## 3. Peptide ProteinAtlas expression values 
## 4. Peptide-MHC binding analysis
## 5. Process public experimental HLA ligand data 
## 6. Analysis and plotting of the results
##
## -----------------------------------------------------------


# 1. Generate recipient-specific peptides from VCF
Rscript ./src/generateRecipPeptides.R "./data/vcf"

# 2.1 Predict immunogenicity of peptides
./src/runImmunogenicity.sh

# 2.2 Select top immunogenic peptides
module load r-env
Rscript ./src/parseImmunogenicity.R

# 3. Map human Protein Atlas expression values to the peptides
module load r-env
Rscript ./src/pepExpression.R

# 4.1 Run peptide binding predictions;
# Parallelized to multiple processors via SLURM batch job system
PATH2='./data/vcf/proteome/immunogenicity'
sbatch ./src/arrayJob_binding.sh

# 4.2 Collect results from peptide info and binding runs
module load r-env
Rscript ./src/collectPepBinding.R

# 4.3 copy result files to ./results subfolder
cp $PATH2/*_peptide.data ./results/peptide_analysis/

# 5.1 Extract and process ligand data
module load r-env
Rscript parse_ligand_data.R
# 5.2 Integrate with recipient-mismatched peptide data
Rscript generate_ligandsets.R

# 6. Regression models and plotting
Rscript ligand_analysis.R





