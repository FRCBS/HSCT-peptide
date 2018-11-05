
## --------------------------------------------------
## Select top immunogenic peptides 
## --------------------------------------------------

library(data.table)
library(Biostrings)

# immunogenicity analysis output folder
setwd('./data/vcf/proteome/immunogenicity')

# path to raw peptide files
pep.path <- './data/vcf/proteome'

# list files
ls.imm <- list.files(pattern='immunogenicity$')
ls.imm.names <- gsub('.9.immunogenicity', '', ls.imm)
ls.pep <- list.files(pep.path, pattern='differential_peptides9.fasta')
ls.pep.names <- gsub('_differential_peptides9.fasta', '', ls.pep)

sapply(1:length(ls.imm), function(i) {
  print(ls.imm.names[i])
  tmp <- fread(ls.imm[i], skip=3, data.table=F) # read immunogenicity analysis results
  tmp <- tmp[tmp[, 3] > 0.2, ] # select top
  ind <- match(ls.imm.names[i], ls.pep.names) # match to original peptides
  tmp2 <- readAAStringSet(paste(pep.path, ls.pep[ind], sep=''))
  pep <- AAStringSet(tmp[, 1])
  tmp2 <- tmp2[match(pep, tmp2)] # same peptides from immunogenicity file
  names(pep) <- paste(names(tmp2), '_', 'Imm=', tmp[, 3], sep='') # ENSEMBL id to top peptides
  Biostrings::writeXStringSet(pep, paste(ls.imm[i], 'top', sep='.')) # write fasta
})

