

library(Biostrings)
library(data.table)
library(stringr)


## Functions

# read binding file, handle missing, sort peptides A-Z
readBinding <- function(x) {
  tt <- tryCatch(fread(x, data.table=F), error=function(x) NA)
  if(!is.na(tt)) tt[order(tt[, 'peptide']), ] else data.frame('percentile_rank'=NA, 'ann_ic50'=NA)
}

# combine peptide information file and binding prediction output files into one
combineFiles <- function(peptide.file, binding.file.A1, binding.file.A2, binding.file.B1, binding.file.B2) {
  data.frame(
    str_split_fixed( names(peptide.file), '_', 11 ),
    rowMeans(cbind(binding.file.A1[, 'percentile_rank'],
                   binding.file.A2[, 'percentile_rank'],
                   binding.file.B1[, 'percentile_rank'],
                   binding.file.B2[, 'percentile_rank']), na.rm=T)
  )
}


# data files folder
setwd('./data/vcf/proteome/immunogenicity')

# Peptide file names
# naming format: GvHD tissues, BM tissues 
# Seq_TranscriptID_GeneID_meanExp_sumExp_ImmunogenicityGvHD_TranscriptID_GeneID_meanExp_sumExp_ImmunogenicityBM
ls.pep <- list.files(pattern='Imm_Exp_GvHDBM.fasta')
ls.pep.nams <- gsub('_Imm_Exp_GvHDBM.fasta', '', ls.pep)

# Peptide binding prediction file names
ls.A1 <- list.files(pattern='A1.IEDB.out')
ls.A1.nams <- gsub('.A1.IEDB.out', '', ls.A1)
ls.A2 <- list.files(pattern='A2.IEDB.out')
ls.A2.nams <- gsub('.A2.IEDB.out', '', ls.A2)
ls.B1 <- list.files(pattern='B1.IEDB.out')
ls.B1.nams <- gsub('.B1.IEDB.out', '', ls.B1)
ls.B2 <- list.files(pattern='B2.IEDB.out')
ls.B2.nams <- gsub('.B2.IEDB.out', '', ls.B2)

# File matching
ls.comm.nams <- Reduce(intersect, list(ls.A1.nams, ls.A2.nams, ls.B1.nams, ls.B2.nams, ls.pep.nams))

# Read in files
pep.files <- lapply(ls.pep[match(ls.comm.nams, ls.pep.nams)], readAAStringSet)
A1.files  <- lapply(ls.A1[match(ls.comm.nams, ls.A1.nams)], readBinding)
A2.files  <- lapply(ls.A2[match(ls.comm.nams, ls.A2.nams)], readBinding)
B1.files  <- lapply(ls.B1[match(ls.comm.nams, ls.B1.nams)], readBinding)
B2.files  <- lapply(ls.B2[match(ls.comm.nams, ls.B2.nams)], readBinding)

# Combine data
combined.files <- lapply(1:length(pep.files), function(i) {
  print(i)
  combineFiles(pep.files[[i]], A1.files[[i]], A2.files[[i]], B1.files[[i]], B2.files[[i]])
})

# write out the data
invisible(sapply(1:length(combined.files), function(i) {
  out.name <- paste(ls.comm.nams[i], 'peptide.data', sep='_')
  write.table(combined.files[[i]], out.name, quote=F, sep='\t', row.names=F,
              col.names=c('pepseq', 'GvHDtranscriptID', 'GvHDgeneID', 'GvHDmeanExp', 'GvHDsumExp', 'GvHDimmunogenicity', 
                          'BMtranscriptID', 'BMgeneID', 'BMmeanExp', 'BMsumExp', 'BMimmunogenicity', 'rank'))
  cat(out.name, sep=' ')
}))


