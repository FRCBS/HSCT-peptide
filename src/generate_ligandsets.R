

## ----------------------------------------------------------------------
##
## Intersect experimentally verified ligand peptides with 
## recipient-mismatched peptide data 
##
## ----------------------------------------------------------------------

library(Biostrings)
library(data.table)
library(tidyverse)


## mHAs

setwd("./data/vcf/proteome")

# read ligand data
mha.lig <- fread('./data/Griffioen_9mers.tsv', header=T, data.table=F)

# read HLA type data
tmp.hla <- read.delim('./data/HLA_types/donor.types', stringsAsFactors=F)[, c('HLA.A', 'HLA.B', 'HLA.C')]
tmp.hla <- apply(tmp.hla, 1, paste, collapse=', ') %>% sapply(., str_split_fixed, ', ',  6) %>% t
tmp.hla <- gsub('HLA-', '', tmp.hla, fixed=T)

# list sample peptide fasta files 
ls.pep <- list.files(pattern='_Imm_Exp_GvHDBM.fasta')
names.pep <- paste0('DT', gsub('_Imm_Exp_GvHDBM.fasta', '', ls.pep, fixed=T))
names.pep <- gsub('DTT', 'DT', names.pep)

# HLA type frequencies
tmp.hla.frq <- table(tmp.hla[tmp.hla!='']) / length(ls.pep)

# match sample order
tmp.hla <- tmp.hla[match(names.pep, rownames(tmp.hla)), ]

# loop over samples
num.peps <- rep(NA, length(ls.pep))
num.mhas <- rep(NA, length(ls.pep))
mha.matches <- imap(ls.pep, function(s, i) {
  print(s)
  # read sample data
  tmp.tab <- as.vector(readAAStringSet(s))
  num.peps[i] <<- length(tmp.tab)
  # subset mHA table by HLA types in the sample 
  tmp.ind <- mha.lig[, 1] %in% tmp.hla[i, ]
  num.mhas[i] <<- sum(tmp.ind)
  # search matches from mHA's
  length(intersect(mha.lig[tmp.ind, 2], tmp.tab))
})

# write output
mha.matches <- data.frame(Sample     = names.pep, 
                          Matches    = mha.matches %>% unlist,
                          PepNum     = num.peps,
                          HLAs       = apply(tmp.hla, 1, function(i) length(unique(i[i!='']))) %>% unname,
                          HLAfreqsum = apply(tmp.hla, 1, function(i) sum(tmp.hla.frq[ names(tmp.hla.frq) %in% unique(i) ])) %>% unname,
                          mHAs       = num.mhas)
write.table(mha.matches, 'mHA_matches.txt', quote=F, sep='\t', row.names=F)





## IEDB ligands

setwd("./data/vcf/proteome")

# read ligand data
iedb.lig <- fread('~/IEDB_ligand.txt', header=F, data.table=F)

# read HLA type data
tmp.hla <- read.delim('./data/HLA_types/donor.types', stringsAsFactors=F)[, c('HLA.A', 'HLA.B', 'HLA.C')]
tmp.hla <- apply(tmp.hla, 1, paste, collapse=', ') %>% sapply(., str_split_fixed, ', ',  6) %>% t
tmp.hla <- gsub('HLA-', '', tmp.hla, fixed=T)

# list sample peptide fasta files 
ls.pep <- list.files(pattern='_differential_peptides9.fasta')
names.pep <- paste0('DT', gsub('_differential_peptides9.fasta', '', ls.pep, fixed=T))
names.pep <- gsub('DTT', 'DT', names.pep)

# HLA type frequencies
tmp.hla.frq <- table(tmp.hla[tmp.hla!='']) / length(ls.pep)

# match sample order
tmp.hla <- tmp.hla[match(names.pep, rownames(tmp.hla)), ]

# loop over samples
num.peps <- rep(NA, length(ls.pep))
num.mhas <- rep(NA, length(ls.pep))
mha.matches <- imap(ls.pep, function(s, i) {
  print(s)
  # read sample data
  tmp.tab <- as.vector(readAAStringSet(s))
  num.peps[i] <<- length(tmp.tab)
  # subset mHA table by HLA types in the sample 
  tmp.ind <- iedb.lig[, 1] %in% tmp.hla[i, ]
  num.mhas[i] <<- sum(tmp.ind)
  # search matches from mHA's
  length(intersect(iedb.lig[tmp.ind, 2], tmp.tab))
})

# write output
mha.matches <- data.frame(Sample     = names.pep, 
                          Matches    = mha.matches %>% unlist,
                          PepNum     = num.peps,
                          HLAs       = apply(tmp.hla, 1, function(i) length(unique(i[i!='']))) %>% unname,
                          HLAfreqsum = apply(tmp.hla, 1, function(i) sum(tmp.hla.frq[ names(tmp.hla.frq) %in% unique(i) ])) %>% unname,
                          mHAs       = num.mhas)
write.table(mha.matches, 'IEDB_Ligand_matches.txt', quote=F, sep='\t', row.names=F)





##  IEDB ligands in immunogenic & expression filtered recipient peptides 

setwd("./data/vcf/proteome/immunogenicity")

# read ligand data
iedb.lig <- fread('~/IEDB_ligand.txt', header=F, data.table=F)

# read HLA type data
tmp.hla <- read.delim('/proj/2000066/McGill_merged/HLA_types/2c.donor.types', stringsAsFactors=F)[, c('HLA.A', 'HLA.B', 'HLA.C')]
tmp.hla <- apply(tmp.hla, 1, paste, collapse=', ') %>% sapply(., str_split_fixed, ', ',  6) %>% t
tmp.hla <- gsub('HLA-', '', tmp.hla, fixed=T)

# list sample peptide fasta files 
ls.pep <- list.files(pattern='_Imm_Exp_GvHDBM.fasta')
names.pep <- paste0('DT', gsub('_Imm_Exp_GvHDBM.fasta', '', ls.pep, fixed=T))
names.pep <- gsub('DTT', 'DT', names.pep)

# HLA type frequencies
tmp.hla.frq <- table(tmp.hla[tmp.hla!='']) / length(ls.pep)

# match sample order
tmp.hla <- tmp.hla[match(names.pep, rownames(tmp.hla)), ]

# loop over samples
num.peps <- rep(NA, length(ls.pep))
num.mhas <- rep(NA, length(ls.pep))
mha.matches <- imap(ls.pep, function(s, i) {
  print(s)
  # read sample data
  tmp.tab <- as.vector(readAAStringSet(s))
  num.peps[i] <<- length(tmp.tab)
  # subset mHA table by HLA types in the sample 
  tmp.ind <- iedb.lig[, 1] %in% tmp.hla[i, ]
  num.mhas[i] <<- sum(tmp.ind)
  # search matches from mHA's
  length(intersect(iedb.lig[tmp.ind, 2], tmp.tab))
})

# write output
mha.matches <- data.frame(Sample     = names.pep, 
                          Matches    = mha.matches %>% unlist,
                          PepNum     = num.peps,
                          HLAs       = apply(tmp.hla, 1, function(i) length(unique(i[i!='']))) %>% unname,
                          HLAfreqsum = apply(tmp.hla, 1, function(i) sum(tmp.hla.frq[ names(tmp.hla.frq) %in% unique(i) ])) %>% unname,
                          mHAs       = num.mhas)
write.table(mha.matches, 'IEDB_Ligand_matches_imm.txt', quote=F, sep='\t', row.names=F)


