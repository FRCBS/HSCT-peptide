

## ---------------------------------------------
## Analysis of ligand count data
## ---------------------------------------------

library(data.table)
library(tidyverse)
library(ggpubr)
library(scales)
library(MASS)

setwd("./data/vcf/proteome")


## Data

iedb.matches.allpep <- fread("./data/IEDB_Ligand_matches.txt", data.table=F))
iedb.matches.allpep <- iedb.matches.allpep[!iedb.matches.allpep$Matches==0, ]
iedb.matches.allpep <- dplyr::rename(iedb.matches.allpep, ID=Sample)

iedb.matches.immpep <- fread("./data/IEDB_Ligand_matches_imm.txt", data.table=F)
iedb.matches.immpep <- dplyr::rename(iedb.matches.immpep, ID=Sample)

mha.matches.allpep  <- fread("./data/IEDB_mHA_matches.txt", data.table=F)
mha.matches.allpep  <- dplyr::rename(mha.matches.allpep, ID=Sample)

aff.ranks.immpep    <- fread("./data/Ranks.tsv", data.table=F)
aff.ranks.immpep$ID <- gsub('T', 'DT', aff.ranks.immpep$ID)


## Parse clinical data

clin    <- read.delim('./data/kaikki_kliiniset_tiedot.txt', stringsAsFactors=F)
clin$ID <- gsub('T', 'DT', clin$ID)
clin    <- clin[match(iedb.matches.allpep$ID, clin$ID), ]
clin    <- clin[!is.na(clin$ID), ]

iedb.matches.allpep <- iedb.matches.allpep[iedb.matches.allpep$ID %in% clin$ID, ]
iedb.matches.immpep <- iedb.matches.immpep[iedb.matches.immpep$ID %in% clin$ID, ]
mha.matches.allpep  <- mha.matches.allpep[mha.matches.allpep$ID   %in% clin$ID, ]
mha.matches.immpep  <- mha.matches.immpep[mha.matches.immpep$ID   %in% clin$ID, ]
aff.ranks.immpep    <- aff.ranks.immpep[aff.ranks.immpep$ID %in% clin$ID, ]

colnames(clin)[11] <- 'relapse'
colnames(clin)[6]  <- 'aGvHD'
colnames(clin)[10] <- 'cGvHD'

metad2 <- read.delim('./data/abc_tiedot.txt', stringsAsFactors=F)
tr.dir <- cbind(metad2[match(iedb.matches.allpep$ID, gsub('SPR.', '', metad2$code)), 'sex'], # Donors
                metad2[match(iedb.matches.allpep$ID %>% gsub('DT', 'T', .), gsub('SPR.', '', metad2$code)), 'sex']) %>% 
  apply(., 1, paste0, collapse='')
tr.dir <- ifelse(tr.dir=='FM', 'FM', 'O') %>% factor         

metad <- read.delim('./data/abc_metadata', stringsAsFactors=F)
metad[, 1] <- gsub('SPR.T', 'DT', metad[, 1])
metad <- metad[match(iedb.matches.allpep$ID, metad[, 1]), ]
clin <- cbind(clin, dage=metad[, 'dage'])

clin[, 'aGvHD'] <- sapply(clin[, 'aGvHD'], function(x) if(x==7 & !is.na(x)) 0 else x)
clin[, 'aGvHD'] <- sapply(clin[, 'aGvHD'], function(x) if(x>=3 & !is.na(x)) '3 or 4' else x)
clin[, 'cGvHD'] <- sapply(clin[, 'cGvHD'], function(x) if(x==7 & !is.na(x)) 0 else x)
clin[which(clin[, 'cGvHD']=='ex ennen'), 10] <- NA
clin$cGvHD <- replace(clin$cGvHD, clin$cGvHD==0, 'No')
clin$cGvHD <- replace(clin$cGvHD, clin$cGvHD==1, 'Limited')
clin$cGvHD <- replace(clin$cGvHD, clin$cGvHD==2, 'Extensive')
clin[which(clin[, 'relapse']=='ex ennen relapse'), 'relapse'] <- NA
clin[, 'relapse'] <- replace(clin$relapse, clin$relapse==1, 'No')
clin[, 'relapse'] <- replace(clin$relapse, clin$relapse==2, 'Yes')

aff.ranks.immpep.ind <- clin$ID %in% aff.ranks.immpep$ID

# create binary phenotypes
acute.bi01 <- sapply(clin$aGvHD, function(x) if((x==3|x==4) & !is.na(x)) 1 else 0)
acute.bi01[is.na(clin$aGvHD)] <- NA
acute.bi02 <- sapply(clin$aGvHD, function(x) if((x==1|x==2|x==3|x==4) & !is.na(x)) 1 else 0)
acute.bi02[is.na(clin$aGvHD)] <- NA
acute.bi03 <- sapply(clin$aGvHD, function(x) if((x==3|x==4) & !is.na(x)) 1 else 0)
acute.bi03[is.na(clin$aGvHD)] <- NA
acute.bi03[(clin$aGvHD==1|clin$aGvHD==2)] <- NA
acute.bi04 <- sapply(clin$aGvHD, function(x) if((x>1) & !is.na(x)) 1 else 0)
acute.bi04[(clin$aGvHD==1)] <- NA
acute.bi04[is.na(clin$aGvHD)] <- NA

chronic.bi01 <- sapply(clin$cGvHD, function(x) if((x==2) & !is.na(x)) 1 else 0)
chronic.bi01[is.na(clin$cGvHD)] <- NA
chronic.bi02 <- sapply(clin$cGvHD, function(x) if((x==1|x==2) & !is.na(x)) 1 else 0)
chronic.bi02[is.na(clin$cGvHD)] <- NA

chronic.bi03 <- clin$cGvHD
chronic.bi03[clin$cGvHD=='Limited']   <- NA
chronic.bi03[clin$cGvHD=='No']   <- 0
chronic.bi03[clin$cGvHD=='Extensive'] <- 1



## Summary stats

tmp.stats <- c(AllPeptides_=iedb.matches.allpep$PepNum,
               ImmPeptides_=iedb.matches.immpep$PepNum,
               AllPeptidesIEDB_=iedb.matches.allpep$Matches,
               ImmPeptidesIEDB_=iedb.matches.immpep$Matches,
               AllPeptidesmHA_=mha.matches.allpep$MHAs,
               HighAffinity_=aff.ranks.immpep$RANK4)
tmp.stats <- data.frame(Measurement=names(tmp.stats) %>% str_split_fixed(., pattern='_', 3) %>% .[, 1] ,
                        Count=tmp.stats, stringsAsFactors=F)

tmp.stats.names <- c(`AllPeptides`    = "All peptides",
                     `AllPeptidesIEDB`= "All peptides & IEDB",
                     `AllPeptidesmHA` = "All peptides & mHAs",
                     `HighAffinity`   = "Predicted high affinity", 
                     `ImmPeptides`    = "Predicted immunogenic peptides",
                     `ImmPeptidesIEDB`= "Predicted immunogenic peptides\n& IEDB")

pdf('./results/Supplementary_Figure_1_count_histograms.pdf', height=7, width=7)
ggplot(tmp.stats, aes(x=Count)) +
  geom_histogram() +
  facet_wrap(tmp.stats$Measurement %>% factor, scales='free', labeller=as_labeller(tmp.stats.names)) + 
  xlab('number of peptides') +
  theme_light() +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1))
dev.off()



## Exploration of ligands counts by clinical outcome

tmp <- full_join(iedb.matches.allpep, clin)[, c('ID', 'Matches', 'aGvHD', 'cGvHD', 'relapse')] %>% 
  dplyr::rename(., All_ligands=Matches) %>% 
  full_join(., iedb.matches.immpep[, c('ID', 'Matches')]) %>% dplyr::rename(., Immunogenic_ligands=Matches) %>% 
  full_join(., mha.matches.allpep[, c('ID', 'MHAs')]) %>% dplyr::rename(., mHAs=MHAs) %>% 
  full_join(., aff.ranks.immpep[, c('ID', 'RANK4')]) %>% dplyr::rename(., Rank=RANK4) %>% 
  dplyr::select(., -ID)

tmp <- gather(tmp, Ligand, Count, c(1, 5:7))

facet.labs <- c(`All_ligands`="M1", # All HLA ligands
                `Immunogenic_ligands`="M2", # Immunogenic HLA ligands
                `mHAs`="M3", # minor H antigens
                `Rank`="M4") # Predicted affinity rank < 4

# aGvHD
p1 <- ggplot(tmp %>% na.omit, aes(x=aGvHD %>% factor, y=Count)) +
  geom_boxplot(outlier.shape=NA) + 
  geom_point(size=0.8, shape=21, color='grey25', alpha=0.5,  position=position_jitter(0.15)) +
  facet_wrap(~Ligand, scales='free', labeller=as_labeller(facet.labs)) +
  theme_light() +
  xlab('aGvHD gradus') +
  ylab('Ligand count') +
  ggtitle('Acute GvHD') +
  theme(strip.background=element_rect(colour="grey70"), plot.title=element_text(size=12))
# cGvHD
p2 <- ggplot(tmp %>% na.omit, aes(x=cGvHD %>% factor(levels=c('No','Limited', 'Extensive'), ordered=T), y=Count)) +
  geom_boxplot(outlier.shape=NA) + 
  geom_point(size=0.8, shape=21, color='grey25', alpha=0.5,  position=position_jitter(0.15)) +
  facet_wrap(~Ligand, scales='free', labeller=as_labeller(facet.labs)) +
  theme_light() +
  xlab('cGvHD gradus') +
  ylab('Ligand count') +
  ggtitle('Chronic GvHD') +
  theme(strip.background=element_rect(colour="grey70"), plot.title=element_text(size=12))
# Relapse
p3 <- ggplot(tmp %>% na.omit, aes(x=relapse %>% factor, y=Count)) +
  geom_boxplot(outlier.shape=NA) + 
  geom_point(size=0.8, shape=21, color='grey25', alpha=0.5,  position=position_jitter(0.15)) +
  facet_wrap(~Ligand, scales='free', labeller=as_labeller(facet.labs)) +
  theme_light() +
  xlab('Relapse occurrence') +
  ylab('Ligand count') +
  ggtitle('Relapse') +
  theme(strip.background=element_rect(colour="grey70"), plot.title=element_text(size=12))

# plot
pdf(height=5, width=10, './results/Figure2_Ligand_counts_by_outcome.pdf')
ggarrange(p1, p3, align='v', nrow=1, ncol=2, labels=c('a', 'b'), font.label=list(size=16))
dev.off()



## Regression models


# chronic GvHD
all.cGvHD.glm <- glm(cGvHD ~ Matches+PepNum+HLAfreqsum+Dir+Age+HLAs, data=data.frame(cGvHD=chronic.bi03, 
                                                                                     iedb.matches.allpep %>% dplyr::select(., -ID) %>% scale,
                                                                                     Age=clin$dage %>% scale, 
                                                                                     Dir=tr.dir), family='binomial') %>% summary
imm.cGvHD.glm <- glm(cGvHD ~ Matches+PepNum+HLAfreqsum+Dir+Age+HLAs, data=data.frame(cGvHD=chronic.bi03, 
                                                                                     iedb.matches.immpep %>% dplyr::select(., -ID) %>% scale,
                                                                                     Age=clin$dage %>% scale, 
                                                                                     Dir=tr.dir), family='binomial') %>% summary
mha.cGvHD.glm <- glm(cGvHD ~ MHAs+PepNum+HLAfreqsum+Dir+Age+HLAs, data=data.frame(cGvHD=chronic.bi03, 
                                                                                  mha.matches.allpep %>% dplyr::select(., -ID) %>% scale,
                                                                                  Age=clin$dage %>% scale, 
                                                                                  Dir=tr.dir), family='binomial') %>% summary
aff.cGvHD.glm <- glm(cGvHD ~ RANK+PepNum+HLAfreqsum+Dir+Age+HLAs, data=data.frame(cGvHD=chronic.bi03[aff.ranks.immpep.ind], 
                                                                   iedb.matches.allpep[aff.ranks.immpep.ind, 3:5] %>% scale,
                                                                   RANK=aff.ranks.immpep$RANK4 %>% scale, 
                                                                   Age=clin$dage[aff.ranks.immpep.ind] %>% scale, 
                                                                   Dir=tr.dir[aff.ranks.immpep.ind]), family=binomial) %>% summary


# acute GvHD
all.aGvHD.polr <- polr(factor(aGvHD) ~ Matches+PepNum+HLAfreqsum+Dir+Age+HLAs, 
                       data=data.frame(aGvHD=clin$aGvHD, 
                                       iedb.matches.allpep %>% dplyr::select(., -ID) %>% scale,
                                       Age=clin$dage %>% scale, 
                                       Dir=tr.dir), Hess=T) %>% summary
all.aGvHD.polr.perm <- replicate(1000, 
                                 summary(polr(factor(aGvHD) ~ Matches+PepNum+HLAfreqsum+Dir+Age+HLAs, 
                                              data=data.frame(aGvHD=sample(clin$aGvHD, 157, F), 
                                                              iedb.matches.allpep %>% dplyr::select(., -ID) %>% scale,
                                                              Age=clin$dage %>% scale, 
                                                              Dir=tr.dir), Hess=T))$coefficients)
all.aGvHD.polr <- data.frame(all.aGvHD.polr$coefficients[1:6, ], perm.p.val=sapply(1:6, function(i) {
  sum(abs(all.aGvHD.polr.perm[i, 3, ]) > abs(all.aGvHD.polr$coefficients[i, 3]))/1000
}))

imm.aGvHD.polr <- polr(factor(aGvHD) ~ Matches+PepNum+HLAfreqsum+Dir+Age+HLAs, 
                       data=data.frame(aGvHD=clin$aGvHD, 
                                       iedb.matches.immpep %>% dplyr::select(., -ID) %>% scale,
                                       Age=clin$dage %>% scale, 
                                       Dir=tr.dir), Hess=T) %>% summary
imm.aGvHD.polr.perm <- replicate(1000, 
                                 summary(polr(factor(aGvHD) ~ Matches+PepNum+HLAfreqsum+Dir+Age+HLAs, 
                                              data=data.frame(aGvHD=sample(clin$aGvHD, 157, F), 
                                                              iedb.matches.immpep %>% dplyr::select(., -ID) %>% scale,
                                                              Age=clin$dage %>% scale, 
                                                              Dir=tr.dir), Hess=T))$coefficients)
imm.aGvHD.polr <- data.frame(imm.aGvHD.polr$coefficients[1:6, ], perm.p.val=sapply(1:6, function(i) {
  sum(abs(imm.aGvHD.polr.perm[i, 3, ]) > abs(imm.aGvHD.polr$coefficients[i, 3]))/1000
}))

mha.aGvHD.polr <- polr(factor(aGvHD) ~ MHAs+PepNum+HLAfreqsum+Dir+Age+HLAs, 
                       data=data.frame(aGvHD=clin$aGvHD, 
                                       mha.matches.allpep %>% dplyr::select(., -ID) %>% scale,
                                       Age=clin$dage %>% scale, 
                                       Dir=tr.dir), Hess=T) %>% summary
mha.aGvHD.polr.perm <- replicate(1000, 
                                 summary(polr(factor(aGvHD) ~ MHAs+PepNum+HLAfreqsum+Dir+Age+HLAs, 
                                              data=data.frame(aGvHD=sample(clin$aGvHD, 157, F), 
                                                              mha.matches.allpep %>% dplyr::select(., -ID) %>% scale,
                                                              Age=clin$dage %>% scale, 
                                                              Dir=tr.dir), Hess=T))$coefficients)
mha.aGvHD.polr <- data.frame(mha.aGvHD.polr$coefficients[1:6, ], perm.p.val=sapply(1:6, function(i) {
  sum(abs(mha.aGvHD.polr.perm[i, 3, ]) > abs(mha.aGvHD.polr$coefficients[i, 3]))/1000
}))

aff.aGvHD.polr <- polr(factor(aGvHD) ~ RANK+PepNum+HLAfreqsum+Dir+Age+HLAs, data=data.frame(aGvHD=clin$aGvHD[aff.ranks.immpep.ind], 
                                                                   iedb.matches.allpep[aff.ranks.immpep.ind, 3:5] %>% scale,
                                                                   RANK=aff.ranks.immpep$RANK4 %>% scale, 
                                                                   Age=clin$dage[aff.ranks.immpep.ind] %>% scale, 
                                                                   Dir=tr.dir[aff.ranks.immpep.ind]), Hess=T) %>% summary
aff.aGvHD.polr.perm <- replicate(1000, 
          summary(polr(factor(aGvHD) ~ RANK+PepNum+HLAfreqsum+Dir+Age+HLAs, data=data.frame(aGvHD=sample(clin$aGvHD[aff.ranks.immpep.ind], 155, F), 
                                                                          iedb.matches.allpep[aff.ranks.immpep.ind, 3:5] %>% scale,
                                                                          RANK=aff.ranks.immpep$RANK4 %>% scale, 
                                                                          Age=clin$dage[aff.ranks.immpep.ind] %>% scale, 
                                                                          Dir=tr.dir[aff.ranks.immpep.ind]), Hess=T))$coefficients)
aff.aGvHD.polr <- data.frame(aff.aGvHD.polr$coefficients[1:6, ], perm.p.val=sapply(1:6, function(i) {
  sum(abs(aff.aGvHD.polr.perm[i, 3, ]) > abs(aff.aGvHD.polr$coefficients[i, 3]))/1000
}))


# relapse

all.relapse.glm <- glm(Relapse ~ Matches+PepNum+HLAfreqsum+Dir+Age+HLAs, data=data.frame(Relapse=relapse.bi, 
                                                                                     iedb.matches.allpep %>% dplyr::select(., -ID) %>% scale,
                                                                                     Age=clin$dage %>% scale, 
                                                                                     Dir=tr.dir), family='binomial') %>% summary
imm.relapse.glm <- glm(Relapse ~ Matches+PepNum+HLAfreqsum+Dir+Age+HLAs, data=data.frame(Relapse=relapse.bi, 
                                                                                     iedb.matches.immpep %>% dplyr::select(., -ID) %>% scale,
                                                                                     Age=clin$dage %>% scale, 
                                                                                     Dir=tr.dir), family='binomial') %>% summary
mha.relapse.glm <- glm(Relapse ~ MHAs+PepNum+HLAfreqsum+Dir+Age+HLAs, data=data.frame(Relapse=relapse.bi, 
                                                                                  mha.matches.allpep %>% dplyr::select(., -ID) %>% scale,
                                                                                  Age=clin$dage %>% scale, 
                                                                                  Dir=tr.dir), family='binomial') %>% summary
aff.relapse.glm <- glm(Relapse ~ RANK+PepNum+HLAfreqsum+Dir+Age+HLAs, data=data.frame(Relapse=relapse.bi[aff.ranks.immpep.ind], 
                                                                                  iedb.matches.allpep[aff.ranks.immpep.ind, 3:5] %>% scale,
                                                                                  RANK=aff.ranks.immpep$RANK4 %>% scale, 
                                                                                  Age=clin$dage[aff.ranks.immpep.ind] %>% scale, 
                                                                                  Dir=tr.dir[aff.ranks.immpep.ind]), family=binomial) %>% summary

# test transplantation date, cGvHD
glm(cGvHD~Year+Dir, data=data.frame(Year=str_split_fixed(clin$Txdate, '\\.', 3)[, 3] %>% 
                                  str_split_fixed(., '\\(', 2) %>% .[, 1] %>% 
                                  str_split_fixed(., '\\.', 5) %>% .[, 1] %>% as.numeric,
                                cGvHD=chronic.bi03, Dir=tr.dir),
    family='binomial') %>% summary




## Collect cGvHD glm results

tmp4 <- data.frame(Ligand=c('M1', 'M2', 'M3'),
                   rbind(all.cGvHD.glm$coef[2, c(1, 2, 2, 4)] %>% exp,
                         imm.cGvHD.glm$coef[2, c(1, 2, 2, 4)] %>% exp,
                         mha.cGvHD.glm$coef[2, c(1, 2, 2, 4)] %>% exp))
tmp4[, 4] <- tmp4[, 4]*-1 # error bar down
colnames(tmp4)[2:5] <- c('Coefficient', 'Errp', 'Errm', 'Pvalue')
tmp4$Pvalue <- c(all.cGvHD.glm$coef[2, 4], imm.cGvHD.glm$coef[2, 4], mha.cGvHD.glm$coef[2, 4])

# plot of reg. coefficients
p6 <- ggplot(tmp4, aes(x=Ligand, y=Coefficient)) +
  geom_point(shape=17, size=4, color='red', alpha=0.7) +
  geom_errorbar(aes(ymin=Errm+Coefficient, ymax=Errp+Coefficient), width=0.1, alpha=0.6) +
  ylim(0.15, 3.3) + #ylim(0.15, 0.9) +
  annotate('text', x=tmp4$Ligand, y=0.2, label=paste('p =', round(tmp4$Pvalue, 3)), size=3) +
  ggtitle('cGvHD logistic regression - No vs. Extensive') +
  xlab('Method') +
  ylab('Odds ratio') +
  theme_light() +
  theme(plot.title=element_text(size=12))

# plot cGvHD ligand counts and reg coeffs
pdf(height=4.5, width=10, './results/Figure3_reg_cGvHD.pdf')
ggarrange(p2, p6, align='h', labels=c('a', 'b'), font.label=list(size=16))
dev.off()


                       
## Analysis of HLA type representation
                       
rec.types <- read.delim('./data/2abc.recipient.types', stringsAsFactors=F)
rec.types <- rec.types[match(gsub('DT', 'T', clin$ID), row.names(rec.types)), ]
rec.types <- cbind(separate(rec.types, HLA.A, c('A_1', 'A_2'), ','), 
                   separate(rec.types, HLA.B, c('B_1', 'B_2'), ','), 
                   separate(rec.types, HLA.C, c('C_1', 'C_2'), ',') )[, c(1,2,9,10,17,18)]
rec.types <- apply(rec.types, 2, function(x) gsub('HLA-', '', x, fixed=T)) %>% data.frame(stringsAsFactors=F)
rec.types <- apply(rec.types, 2, function(x) gsub(' ', '', x, fixed=T)) %>% data.frame(stringsAsFactors=F)

iedb.hla <- fread('./data/IEDB_ligand.txt', data.table=F, header=F)
iedb.hla <- iedb.hla[nchar(iedb.hla$V1)>6, ]
griffioen.hla <- fread('./data/Griffioen_9mers.tsv', data.table=F)
griffioen.hla <- griffioen.hla[nchar(griffioen.hla$HLA)>6, ]
  
sum(c(rec.types$A_1, rec.types$A_2) %in% iedb.hla$V1)       / 314
sum(c(rec.types$A_1, rec.types$A_2) %in% griffioen.hla$HLA) / 314
sum(c(rec.types$B_1, rec.types$B_2) %in% iedb.hla$V1)       / 314
sum(c(rec.types$B_1, rec.types$B_2) %in% griffioen.hla$HLA) / 314
sum(c(rec.types$C_1, rec.types$C_2) %in% iedb.hla$V1)       / 314
sum(c(rec.types$C_1, rec.types$C_2) %in% griffioen.hla$HLA) / 314


