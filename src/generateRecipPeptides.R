
## -------------------------------------------------------
## Generate peptides from non-synonymous variants
## that are different between donor and recipient
## -------------------------------------------------------

# libs
library(Biostrings)

# This script is run in the ./src dir
# Load  needed functions and 
source('./src/vcf_functions.R')
grch37.transcripts <- filterTranscripts('./data/Homo_sapiens.GRCh37.75.cds.all.fa')

# Argument is the path to output folder
args <- commandArgs(trailingOnly=T)[1]
setwd(args) # change to the folder given in args

# read VCF variant files
ls.vcf <- list.files('../annotated/', 'vcf.gz$')
ls.vcf.pairids <- intersect(gsub('SPR.T|.filt.anno.vcf.gz', '',  ls.vcf[grepl('SPR.T', ls.vcf)]),
                            gsub('SPR.DT|.filt.anno.vcf.gz', '', ls.vcf[!grepl('SPR.T', ls.vcf)]))
ls.vcf.dpind <- t(sapply(ls.vcf.pairids, function(x) which(grepl(x, ls.vcf)))) # donor-patient sample match table

# loop over files..
for(i in 1:nrow(ls.vcf.dpind)) {
  compareProteomes(paste('../annotated/', ls.vcf[ls.vcf.dpind[i, 1]], sep=''), # donor
                   paste('../annotated/', ls.vcf[ls.vcf.dpind[i, 2]], sep='')) # patient
  # generates peptides, compare, select specific to patient
  don.pep <- generatePeptides2(readAAStringSet(paste('DT', ls.vcf.pairids[i], '_unique_proteins.fasta', sep='')))
  pat.pep <- generatePeptides2(readAAStringSet(paste('T',  ls.vcf.pairids[i], '_unique_proteins.fasta', sep='')))
  pat.pep.uniq <- lapply(1:length(pat.pep), function(x) setdiff(pat.pep[[x]], don.pep[[x]]))
  # write peptides specific to patient to fasta files
  sapply(1:length(pat.pep.uniq), function(x) {
    mer <- c(8, 9, 10, 11, 12, 15)
    out.name <- paste('T', ls.vcf.pairids[i], sep='')
    invisible(writeXStringSet(pat.pep.uniq[[x]], paste(out.name, '_', 'differential_peptides',  mer[x], '.fasta', sep='')))
  })
}


