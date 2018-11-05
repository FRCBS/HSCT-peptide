

## libs

library(seqminer)
library(Biostrings)


## -------------------------------------------------------------------------
## Formatting functions 
## -------------------------------------------------------------------------

# identical ID's' for donors and recipients
formatNam <- function(x) gsub('[SPR.T|SPR.DT|T|DT]|-', '', x) 

# extract ID
extractNam <- function(x) sapply(x, function(y) strsplit(y, '_')[[1]][1]) 

# remove NA's
rmNA <- function(d) d[!is.na(d)] 


## -------------------------------------------------------------------------
## Filter transcripts FASTA sequence file 
## Include only sequences starting with ATG, remove =<50 nt sequences
## -------------------------------------------------------------------------

filterTranscripts <- function(file='Homo_sapiens.GRCh37.75.cds.all.fa') {
  trs <- Biostrings::readDNAStringSet(file)
  trs <- trs[which(Biostrings::width(trs) > 50)]
  start.codon <- Biostrings::DNAString('ATG')
  true.start <- Biostrings::subseq(trs, 1, 3)
  true.start <- which(true.start==start.codon)
  trs <- trs[true.start]
  names(trs) <- sapply(names(trs), function(x) strsplit(x, ' ')[[1]][1])
  trs
}


## -------------------------------------------------------------------------
##
## Functions for filtering VCF variant data
##
## -------------------------------------------------------------------------

# Extract info field from a tabix indexed VCF
getVCFInfo <- function(file) { 
  unlist(sapply(c(1:22, 'X', 'Y'), function(x) {
    seqminer::readVCFToListByRange(file, paste(x, '1-250000000', sep=':'), vcfInfo='', annoType='',
                              vcfColumn="INFO", vcfIndv="GT")$INFO
  }))
}

# Select certain variant effects, input is VCF info field
filterEff2 <- function(info) { 
  sapply(info, function(x) {
    effs <- unlist(strsplit(strsplit(x, 'ANN=')[[1]][2], ',')[[1]])
    effs <- effs[grepl('protein_coding', effs)]
    effs <- effs[!grepl('WARNING_TRANSCRIPT_INCOMPLETE|sequence_feature|intron_variant|regulatory_region_variant|UTR_variant', effs)]
    if(sum(grepl('missense|insertion|deletion|frameshift|stop_gained', effs)) > 0) T else F
  })
}

# Extract transcript ENSEMBL ID from a variant effect
extractTranscript <- function(tt) { 
  out <- strsplit(tt, '|', fixed=T)[[1]] 
  out[which(out == 'transcript') + 1] # transcript ENSMBL ID 
}


## -------------------------------------------------------------------------
##
## Functions for extracting coding sequence variants 
## from SNPEff annotated  VCF data
##
## -------------------------------------------------------------------------

# remove empty items from a vector
remEmpty <- function(vec) vec[!(vec == '')] 

# parse missense mutation: |c.442C>T|
parseMis <- function(pat) {
  subst <- remEmpty(strsplit(pat, 'c.', fixed=T)[[1]][2])
  subst <- remEmpty(strsplit(subst, '|', fixed=T)[[1]][1])
  nucl  <- remEmpty(strsplit(subst, '>', fixed=T)[[1]][2]) # nucl from to
  posit <- remEmpty(strsplit(subst, '[A-z]')[[1]][1])
  c(posit, nucl)
}

# parse (inframe) deletion: |c.528_530delCAA| |c.4862delA|
parseDel <- function(pat) {
  prs <- remEmpty(strsplit(pat, 'c.', fixed=T)[[1]][2])
  out <- remEmpty(strsplit(strsplit(prs, 'del', fixed=T)[[1]][1], '_')[[1]]) # return deleted nucl posit
  if(length(out) == 1) c(out, out) else {
    if(is.na(out[1]) | is.na(out[2])) c(NA, NA) else out
  }
}

# parse (in frame) insertion |c.344_349dupAAGAAA| |c.6dupG|
# insert sequence after position
parseIns <- function(pat) {
  subst <- remEmpty(strsplit(pat, 'c.', fixed=T)[[1]][2])
  subst <- remEmpty(strsplit(subst, '|', fixed=T)[[1]][1])
  nucl  <- remEmpty(strsplit(subst, 'dup|ins')[[1]][2]) 
  posit <- as.numeric(remEmpty(strsplit(strsplit(subst, 'dup|ins')[[1]][1], '_')[[1]]))
  posit <- max(posit, na.rm=T)
  c(posit, nucl) # return start position, nucls
}

# parse frameshift insertions |c.6dupG|
parseFrameshiftIns <- function(pat) {
  nucl <- strsplit(pat, 'dup|ins')[[1]][2] # inserted nucl seq
  posit <- strsplit(strsplit(pat, 'ins')[[1]][1], '/c.', fixed=T)[[1]][2] # inserted nucl position
  posit <- strsplit(posit, '_')[[1]][1] # insert new nucleotide after this position
  c(posit, nucl) # return start position, nucls
}

# parse DEL frameshift p.Trp176fs/c.528delG
parseFrameshiftDel <- function(pat) {
  out <- strsplit(strsplit(pat, 'del')[[1]][1], '/c.', fixed=T)[[1]][2] # deleted nucl position(s)
  t(sapply(out, function(x) {
    tt <- strsplit(x, '_')[[1]]
    if(length(tt)==1) c(tt, tt) else tt
  }))
}

# parse gained STOP substitution |c.442C>T|
parseStop <- function(pat) {
  subst <- remEmpty(strsplit(pat, 'c.', fixed=T)[[1]][2])
  subst <- remEmpty(strsplit(subst, '|', fixed=T)[[1]][1])
  nucl  <- remEmpty(strsplit(subst, '>', fixed=T)[[1]][2]) 
  posit <- as.numeric(remEmpty(strsplit(strsplit(subst, 'dup|ins')[[1]][1], '_')[[1]]))
  post  <- max(posit, na.rm=T)
  c(posit, nucl) # return start position, nucls
}

# Extract all coding variants from a filtered VCF
findEffects2 <- function(info) {
  mis.dat <- ins.dat <- del.dat <- fin.dat <- fde.dat <- sto.dat <- c()
  for(i in 1:length(info)) {
    #print(i)
    effs <- unlist(strsplit(strsplit(info[i], 'ANN=')[[1]][2], ',')[[1]])
    effs.types <- unname(unlist(sapply(effs, function(ii) strsplit(ii, '|', fixed=T)[[1]][2])))
    missens <- which(grepl('missense',    effs.types))
    inserti <- which(grepl('insertion',   effs.types))
    deletio <- which(grepl('deletion',    effs.types))
    framein <- which(grepl('frameshift',  effs.types) & grepl('ins|dup',  effs))
    framede <- which(grepl('frameshift',  effs.types) & grepl('del',      effs))
    stopgai <- which(grepl('stop_gained', effs.types))
    if(length(missens)>0) mis.dat <- rbind(mis.dat, t(sapply(missens, function(x) c(parseMis(effs[x]), extractTranscript(effs[x])))))
    if(length(inserti)>0) ins.dat <- rbind(ins.dat, t(sapply(inserti, function(x) c(parseIns(effs[x]), extractTranscript(effs[x])))))
    if(length(deletio)>0) del.dat <- rbind(del.dat, t(sapply(deletio, function(x) c(parseDel(effs[x]), extractTranscript(effs[x])))))
    if(length(framein)>0) fin.dat <- rbind(fin.dat, t(sapply(framein, function(x) c(parseIns(effs[x]), extractTranscript(effs[x])))))
    if(length(framede)>0) fde.dat <- rbind(fde.dat, t(sapply(framede, function(x) c(parseDel(effs[x]), extractTranscript(effs[x])))))
    if(length(stopgai)>0) sto.dat <- rbind(sto.dat, t(sapply(stopgai, function(x) c(parseMis(effs[x]), extractTranscript(effs[x])))))
  }
  list(mis.dat, ins.dat, del.dat, fin.dat, fde.dat, sto.dat)
}


## -------------------------------------------------------------------------
##
## Functions for transferring variant data to transcript data
##
## -------------------------------------------------------------------------

# insertions
doIns <- function(ins, pos, target) {
  if(is.na(ins)) target else {
    ins <- DNAString(ins)
    if(pos < length(target)) {
      c(target[1:pos], ins, target[(as.numeric(as.vector(pos))+1):length(target)]) 
    } else c(target, ins)
  }
}

# substitutions
doSubs <- function(subs, pos, target) {
  replaceLetterAt(target, as.numeric(pos), as.character(subs))
}

# deletions
doDel <- function(from, to, target) {
  if(length(target) <= to) to <- (length(target)-1)
  if(from < 2) {
    target[(to+1):length(target)]
    } else c(target[1:(from-1)], target[(to+1):length(target)]) 
}

# incorporate missense mutations into transcripts
implementMissense <- function(tbl=vcf.eff, transcripts=trs, transcripts.names=trs.names) {
  inds <- match(tbl[[1]][, 3], transcripts.names) # match transcript ids
  tb <- data.frame(tbl[[1]][, 1:3], inds) # merged variant table and transcipt indices
  # tb columns: position, base, id, index
  tb <- tb[!is.na(tb[, 4]), ] # remove if transcript id is not found
  tb <- tb[!duplicated(tb), ]
  tb[, 2] <- as.character(tb[, 2]) # base to char
  tb[, 1] <- as.numeric(as.vector(tb[, 1])) # position to num
  tb <- tb[!is.na(tb[, 1]), ] 
  tb <- tb[!is.na(tb[, 2]), ] 
  # loop over variants
  out <- DNAStringSet(sapply(1:nrow(tb), function(x) {
    doSubs(tb[x, 2], tb[x, 1], transcripts[[tb[x, 4]]])
  }))
  names(out) <- tb[, 3] # IDs as names
  unique(out)
}

# incorporate inframe insertions to transcripts
implementIns <- function(tbl=vcf.eff, transcripts=trs, transcripts.names=trs.names) {
  inds <- match(tbl[[2]][, 3], transcripts.names) # match transcript ids
  tb <- data.frame(tbl[[2]][, 1:3], inds) # merged variant table and transcipt indices
  # tb columns: position, base, id, index
  tb <- tb[!is.na(tb[, 4]), ] # remove if transcript id is not found
  tb <- tb[!duplicated(tb), ]
  tb[, 2] <- as.character(tb[, 2]) # base to char
  tb[, 1] <- as.numeric(as.vector(tb[, 1])) # position to num
  tb <- tb[!is.na(tb[, 1]), ] 
  tb <- tb[!is.na(tb[, 2]), ] 
  # loop over variants
  out <- DNAStringSet(sapply(1:nrow(tb), function(x) {
    doIns(tb[x, 2], tb[x, 1], transcripts[[tb[x, 4]]])
  }))
  names(out) <- tb[, 3] # IDs as names
  unique(out)
}

# incorporate inframe deletions to transcripts
implementDel <- function(tbl=vcf.eff, transcripts=trs, transcripts.names=trs.names) {
  inds <- match(tbl[[3]][, 3], transcripts.names) # match transcript ids
  tb <- data.frame(tbl[[3]][, 1:3], inds) # merged variant table and transcipt indices
  # tb columns: position, base, id, index
  tb <- tb[!is.na(tb[, 4]), ] # remove if transcript id is not found
  tb <- tb[!duplicated(tb), ]
  tb[, 1] <- as.numeric(as.vector(tb[, 1])) # position 'from' to num
  tb[, 2] <- as.numeric(as.vector(tb[, 2])) # position 'to' to num
  tb <- tb[!is.na(tb[, 1]), ] 
  tb <- tb[!is.na(tb[, 2]), ] 
  # loop over variants
  out <- DNAStringSet(sapply(1:nrow(tb), function(x) {
    doDel(tb[x, 1], tb[x, 2], transcripts[[tb[x, 4]]])
  }))
  names(out) <- tb[, 3] # IDs as names
  unique(out)
}

# incorporate frameshift insertions to transcripts
implementFSIns <- function(tbl, transcripts, transcripts.names) {
  inds <- match(tbl[[4]][, 3], transcripts.names) # match transcript ids
  tb <- data.frame(tbl[[4]][, 1:3], inds) # merged variant table and transcipt indices
  # tb columns: position, base, id, index
  tb <- tb[!is.na(tb[, 4]), ] # remove if transcript id is not found
  tb <- tb[!is.na(tb[, 1]), ]
  tb <- tb[!is.na(tb[, 2]), ] 
  tb <- tb[!duplicated(tb), ]
  tb[, 2] <- as.character(tb[, 2]) # base to char
  tb[, 1] <- as.numeric(as.vector(tb[, 1])) # position to num
  tb <- tb[!(tb[, 1] < 1), ] # remove negative positions
  # loop over variants
  out <- DNAStringSet(sapply(1:nrow(tb), function(x) {
    doIns(tb[x, 2], tb[x, 1], transcripts[[tb[x, 4]]])
  }))
  names(out) <- tb[, 3] # IDs as names
  unique(out)
}

# incorporate frameshift deletions to transcripts
implementFSDel <- function(tbl=vcf.eff, transcripts=trs, transcripts.names=trs.names) {
  inds <- match(tbl[[5]][, 3], transcripts.names) # match transcript ids
  tb <- data.frame(tbl[[5]][, 1:3], inds) # merged variant table and transcipt indices
  # tb columns: position, base, id, index
  tb <- tb[!is.na(tb[, 4]), ] # remove if transcript id is not found
  tb <- tb[!duplicated(tb), ]
  tb[, 1] <- as.numeric(as.vector(tb[, 1])) # position 'from' to num
  tb[, 2] <- as.numeric(as.vector(tb[, 2])) # position 'to' to num
  tb <- tb[!is.na(tb[, 1]), ] 
  tb <- tb[!is.na(tb[, 2]), ] 
  # loop over variants
  out <- DNAStringSet(sapply(1:nrow(tb), function(x) {
    doDel(tb[x, 1], tb[x, 2], transcripts[[tb[x, 4]]])
  }))
  names(out) <- tb[, 3] # IDs as names
  unique(out)
}

# incorporate gained stops to transcripts
implementGStop <- function(tbl=vcf.eff, transcripts=trs, transcripts.names=trs.names) {
  inds <- match(tbl[[6]][, 3], transcripts.names) # match transcript ids
  tb <- data.frame(tbl[[6]][, 1:3], inds) # merged variant table and transcipt indices
  # tb columns: position, base, id, index
  tb <- tb[!is.na(tb[, 4]), ] # remove if transcript id is not found
  tb <- tb[!duplicated(tb), ]
  tb[, 2] <- as.character(tb[, 2]) # base to char
  tb[, 1] <- as.numeric(as.vector(tb[, 1])) # position to num
  tb <- tb[!is.na(tb[, 1]), ] 
  tb <- tb[!is.na(tb[, 2]), ] # remove if NA base
  # loop over variants
  out <- DNAStringSet(sapply(1:nrow(tb), function(x) {
    doSubs(tb[x, 2], tb[x, 1], transcripts[[tb[x, 4]]])
  }))
  names(out) <- tb[, 3] # IDs as names
  unique(out)
}

# output all variant transcripts as a stringset
implementVariants <- function(variant.table, transcript.data) {
  
  missens.seq <- implementMissense(variant.table, transcripts=transcript.data, transcripts.names=names(transcript.data))
  ins.seq     <- implementIns(variant.table,      transcripts=transcript.data, transcripts.names=names(transcript.data))
  del.seq     <- implementDel(variant.table,      transcripts=transcript.data, transcripts.names=names(transcript.data))
  fsins.seq   <- implementFSIns(variant.table,    transcripts=transcript.data, transcripts.names=names(transcript.data))
  fsdel.seq   <- implementFSDel(variant.table,    transcripts=transcript.data, transcripts.names=names(transcript.data))
  gstop.seq   <- implementGStop(variant.table,    transcripts=transcript.data, transcripts.names=names(transcript.data))
  
  out <- list(missens.seq, ins.seq, del.seq, fsins.seq, fsdel.seq, gstop.seq)
  out.names <- unlist(sapply(out, names))
  out <- do.call(c, out)
  names(out) <- out.names
  out
}



## -------------------------------------------------------------------------
##
## Functions for generating HSCT recipient-mismatched peptides 
##
## -------------------------------------------------------------------------

# produce all n-mer peptides (len) from a protein sequence (prot; AAString)
toPep <- function(prot, len) {
  AAStringSet(Views(prot, 1:(length(prot) - (len-1)), len:length(prot)))
}

# generate n-mer peptides for each protein from AAStringSet; 
# preserve transcript IDs and generate peptide names
proteomeToPep <- function(proteome, len) {
  out <- lapply(proteome, toPep, len=len) # list of AAStringSets for all peptides per protein
  out.lengths <- sapply(out, length) # number of peptides in each element of 'out'
  out.names <- rep(names(out.lengths), times=out.lengths)
  out.names <- paste(out.names, unlist(sapply(out.lengths, function(x) paste('pep', 1:x, sep=''))), sep='_') # peptide nums
  out <- do.call(c, unname(out)) # collapse the list
  names(out) <- out.names # name the AAStringSet
  out
}

# translate into protein
translateRNA <- function(nucl) {
  non.atg.ind <- which(subseq(nucl, 1, 3) != DNAString('ATG'))
  if(length(non.atg.ind) > 0) {
    for(pp in non.atg.ind) { # find new start
      non.atg   <- nucl[[pp]]
      new.start <- which(DNAStringSet(sapply(1:(length(non.atg)-3), function(x) subseq(non.atg, x, x+2))) == DNAStringSet('ATG'))[1]
      if(!is.na(new.start)) { # if found, start sequence from that
        nucl[[pp]]  <- non.atg[new.start : length(non.atg)]
      } else nucl[[pp]] <- DNAString('A') # replace NA's with single A
    }
  }
  nucl <- nucl[!(nucl == DNAString('A'))] # remove NA's
  prot <- suppressWarnings(translate(nucl)) # translation
  prot <- AAStringSet(gsub('*', '', prot, fixed=T)) # protein sequences without stop symbol (*)
  prot 
}

# compare proteomes between donor and recipient
compareProteomes <- function(vcf.donor, vcf.patient, ref.transcripts=grch37.transcripts) {
  # output ID's
  don.name <- strsplit(vcf.donor, '.', fixed=T)[[1]][4]
  pat.name <- strsplit(vcf.patient, '.', fixed=T)[[1]][4]
  # generate donor proteome
  don <- getVCFInfo(vcf.donor)
  don <- don[grepl('missense|insertion|deletion|frameshift|stop_gained', don)]
  don <- don[which(filterEff2(don))]
  don <- suppressWarnings(findEffects2(don))
  don <- suppressWarnings(implementVariants(don, ref.transcripts))
  don <- translateRNA(don)
  # generate patient proteome
  pat <- getVCFInfo(vcf.patient)
  pat <- pat[grepl('missense|insertion|deletion|frameshift|stop_gained', pat)]
  pat <- pat[which(filterEff2(pat))]
  pat <- suppressWarnings(findEffects2(pat))
  pat <- suppressWarnings(implementVariants(pat, ref.transcripts))
  pat <- translateRNA(pat)
  # identical proteins between patient and donor
  iden <- intersect(pat, don) 
  # remove identical proteins 
  don <- don[-na.omit(match(names(iden), names(don)))]
  pat <- pat[-na.omit(match(names(iden), names(pat)))]
  # write fasta
  writeXStringSet(don, paste(don.name, 'unique_proteins.fasta', sep='_'))
  writeXStringSet(pat, paste(pat.name, 'unique_proteins.fasta', sep='_'))
}


# generate peptides using Biostrings functions
# input is a DNAStringSet containing the transcriptome
generatePeptides2 <- function(input, mer=c(8, 9, 10, 11, 12, 15), out.name=NULL, writefasta=F, translate=F) {
  # translate transcripts
  if(translate) prot <- translateRNA(input) else prot <- input
  # generate peptides of different lengths given by 'mer'
  out <- lapply(mer, function(x) proteomeToPep(prot, len=x))
  # write peptide fasta files
  if(writefasta & !is.null(out.name)) {
    sapply(1:length(out), function(x) writeXStringSet(out[[x]], paste(out.name, '_', mer[x], 'peptides.fasta', sep='')))
  } else out
}  



