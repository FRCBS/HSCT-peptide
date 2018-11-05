
## libs
library(data.table) 
library(dplyr) 
library(Biostrings)


setwd('./data')


## ---------------------------------------------
## Parse mHA peptide data
## ---------------------------------------------

tmp <- fread('./data/Griffioen_mHAs.tsv', data.table=F)
tmp <- data.frame(Peptide=c(sapply(tmp$mHA, function(i) gsub('\\[./', '', i) %>% gsub('\\]', '', .)), 
                            sapply(tmp$mHA, function(i) gsub('/.\\]', '', i) %>% gsub('\\[', '', .))),
                  HLA=c(tmp$HLA, tmp$HLA), Gene=c(tmp$Gene, tmp$Gene), 
                  stringsAsFactors=F)
tmp2 <- sapply(tmp$Peptide, function(i) str_sub(i, 1:3, 9:11)) %>% t
tmp <- data.frame(Peptide=c(tmp2[, 1], tmp2[, 2], tmp2[, 3]), 
                  HLA=c(tmp$HLA, tmp$HLA, tmp$HLA), Gene=c(tmp$Gene, tmp$Gene, tmp$Gene), 
                  stringsAsFactors=F)
tmp <- tmp[!sapply(tmp$Peptide, nchar) != 9, ]
tmp <- tmp[!tmp %>% duplicated, ]
write.table(tmp, './data/Griffioen_9mers.tsv', sep='\t', quote=F, row.names=F)


## ---------------------------------------------
## Parse IEDB HLA ligand data
## ---------------------------------------------

tmp <- fread('mhc_ligand_full.csv', data.table=F, skip=1)
tmp <- tmp[which(tmp[, 'Organism Name']=='Homo sapiens' & 
                   (tmp[, 'Qualitative Measure']=='Positive' | tmp[, 'Qualitative Measure']=='Positive-High') &
                   tmp[, 'Name']=='Homo sapiens'), ]
tmp <- tmp[grepl('HLA', tmp[, 'Allele Name']), ]
tmp <- tmp[!grepl('class', tmp[, 'Allele Name']), ]
tmp <- tmp[!grepl('mutant', tmp[, 'Allele Name']), ]
tmp <- tmp[, c('Allele Name', 'Description')]
tmp <- tmp[!duplicated(tmp), ]
tmp[, 1] <- gsub('HLA-', '', tmp[, 1], fixed=T)
colnames(tmp) <- c('Allele', 'Peptide')
tmp[, 2] <- sapply(tmp[, 2], function(x) strsplit(x, ' ')[[1]][1])
write.table(tmp, 'IEDB_ligand.txt', quote=F, sep='\t', row.names=F, col.names=F)



