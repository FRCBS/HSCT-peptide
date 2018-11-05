## -----------------------------------------------------
## Mapping expression data to peptides
## -----------------------------------------------------

library(Biostrings)
library(data.table)
library(stringr)
options(stringsAsFactors=F)

# Human Protein Atlas tissue expression data
normal.tissue <- read.csv('~/normal_tissue.csv', stringsAsFactors=F) # The Human Protein Atlas version 16 and Ensembl version 83.38
normal.tissue <- normal.tissue[normal.tissue[, 'Reliability']=='Supportive', ]
normal.tissue <- normal.tissue[!normal.tissue[, 'Level']=='Not detected', ]
normal.tissue.gvhd <- normal.tissue[normal.tissue[, 'Tissue'] %in% c('colon', 'duodenum', 'small intestine', 'liver', 'lung', 'skin 1', 'skin 2'), ]
normal.tissue.bm   <- normal.tissue[normal.tissue[, 'Tissue'] %in% 'bone marrow', ]
bm.geneids <- intersect(names((which(table(normal.tissue[, 1])==1))), normal.tissue.bm[, 1]) # strictly bm specific genes
ensembl <- read.delim('~/ENSEMBL_IDs.txt', header=T, stringsAsFactors=F) 


# Immunogenicity analysis output folder
setwd('./data/vcf/proteome/immunogenicity')

# Peptide file list
ls.imm       <- list.files(pattern='9.immunogenicity.top') # immunogenicity filtered peptide files
ls.imm.names <- gsub('.9.immunogenicity.top', '', ls.imm) # sample names

# Map expression information to peptides
sapply(1:length(ls.imm), function(i) {

  print(ls.imm.names[i])

  tmp <- readAAStringSet(ls.imm[i]) # read peptide file
  trs <- str_split(names(tmp), '_', n=2, T)[, 1]  # extract ensembl transcript id
  ids <- ensembl[match(trs, ensembl[, 2]), 1] # map to gene ids

  # GvHD tissue expression
  normal.tissue.gvhd.sel <- normal.tissue.gvhd[normal.tissue.gvhd[, 1] %in% ids, ] # select from protein expression data
  normal.tissue.gvhd.sel[, 'Level'] <- sapply(normal.tissue.gvhd.sel[, 'Level'], function(x) switch(x, Low=1, Medium=2, High=3)) # to numerical
  gene.mean <- tapply(normal.tissue.gvhd.sel[, 'Level'], normal.tissue.gvhd.sel[, 'Gene'], mean) # mean expr level per gene
  gene.sum  <- tapply(normal.tissue.gvhd.sel[, 'Level'], normal.tissue.gvhd.sel[, 'Gene'], sum) # sum expr level per gene
  ind <- merge(data.frame(Gene=names(gene.sum), ExpSum=gene.sum), ensembl, by.x='Gene', by.y='Ensembl.Gene.ID') # join with ENSEMBL
  ind <- merge(data.frame(Gene=names(gene.mean), ExpMean=gene.mean), ind, by.x='Gene', by.y='Gene')
  ind <- ind[na.omit(match(trs, ind[, 'Ensembl.Transcript.ID'])), ] 
  ind <- ind[!duplicated(ind), ]
  trs.ind <- trs %in% ind[, 'Ensembl.Transcript.ID']
  trs.tmp <- trs[trs.ind] # select transcripts that pass expression filtering
  tmp.gvhd <- tmp[trs.ind] # select transcripts that pass expression filtering
  trs.ind <- merge(data.frame(TRS=trs.tmp), ind, by.x='TRS', by.y='Ensembl.Transcript.ID') # align expression values with all peptides
  trs.ind.gvhd <- trs.ind[match(trs.tmp, trs.ind[, 'TRS']), ]
  trs.ind.gvhd <- data.frame(trs.ind.gvhd, Imm=str_split_fixed(str_split_fixed(names(tmp.gvhd), 'Imm=', n=2)[, 2], '_', n=2)[, 1]) #join immunogenicity

  # Bone marrow tissue expression
  normal.tissue.bm.sel <- normal.tissue.bm[normal.tissue.bm[, 1] %in% ids, ] # select from protein expression data
  normal.tissue.bm.sel[, 'Level'] <- sapply(normal.tissue.bm.sel[, 'Level'], function(x) switch(x, Low=1, Medium=2, High=3)) # to numerical
  gene.mean <- tapply(normal.tissue.bm.sel[, 'Level'], normal.tissue.bm.sel[, 'Gene'], mean) # mean expr level per gene
  gene.sum  <- tapply(normal.tissue.bm.sel[, 'Level'], normal.tissue.bm.sel[, 'Gene'], sum) # sum expr level per gene
  ind <- merge(data.frame(Gene=names(gene.sum), ExpSum=gene.sum), ensembl, by.x='Gene', by.y='Ensembl.Gene.ID') # join with ENSEMBL
  ind <- merge(data.frame(Gene=names(gene.mean), ExpMean=gene.mean), ind, by.x='Gene', by.y='Gene')
  ind <- ind[na.omit(match(trs, ind[, 'Ensembl.Transcript.ID'])), ] 
  ind <- ind[!duplicated(ind), ]
  trs.ind <- trs %in% ind[, 'Ensembl.Transcript.ID']
  trs.tmp <- trs[trs.ind] # select transcripts that pass expression filtering
  tmp.bm <- tmp[trs.ind] # select transcripts that pass expression filtering
  trs.ind <- merge(data.frame(TRS=trs.tmp), ind, by.x='TRS', by.y='Ensembl.Transcript.ID') # align expression values with all peptides
  trs.ind.bm <- trs.ind[match(trs.tmp, trs.ind[, 'TRS']), ]
  trs.ind.bm <- data.frame(trs.ind.bm, Imm=str_split_fixed(str_split_fixed(names(tmp.bm), 'Imm=', n=2)[, 2], '_', n=2)[, 1]) #join immunogenicity

  # merge GvHD and BM peptides
  test.m <- merge(data.frame(tmp.gvhd, trs.ind.gvhd), data.frame(tmp.bm, trs.ind.bm), by.x=1, by.y=1, all=T)
  all.pep <- AAStringSet(test.m[, 1])
  names(all.pep) <- apply(test.m, 1, paste, collapse='_')
  writeXStringSet(all.pep, paste(ls.imm.names[i], '_Imm_Exp_GvHDBM.fasta', sep='')) # write output fasta

})


