

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

iedb.matches.allpep <- fread("./data/IEDB_Ligand_matches.txt", data.table=F)
iedb.matches.allpep <- iedb.matches.allpep[!iedb.matches.allpep$Matches==0, ]
iedb.matches.allpep <- dplyr::rename(iedb.matches.allpep, ID=Sample)

iedb.matches.immpep <- fread("./data/IEDB_Ligand_matches_imm.txt", data.table=F)
iedb.matches.immpep <- dplyr::rename(iedb.matches.immpep, ID=Sample)

mha.matches.allpep  <- fread("./data/IEDB_mHA_matches.txt", data.table=F)
mha.matches.allpep  <- dplyr::rename(mha.matches.allpep, ID=Sample)

aff.ranks.immpep    <- fread("./data/Ranks.tsv", data.table=F)
aff.ranks.immpep$ID <- gsub('T', 'DT', aff.ranks.immpep$ID)


## Parse clinical data

clin    <- read.delim('./data/kaikki_kliiniset_tiedot_1216.txt', stringsAsFactors=F)
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
clinical.cgvhd <- as.numeric(clin$cGvHD %>% factor(., levels=c('No', 'Limited', 'Extensive')))
clinical.cgvhd <- clinical.cgvhd-1

chronic.bi01 <- sapply(clinical.cgvhd, function(x) if((x==2) & !is.na(x)) 1 else 0)
chronic.bi01[is.na(clinical.cgvhd)] <- NA
chronic.bi02 <- sapply(clinical.cgvhd, function(x) if((x==1|x==2) & !is.na(x)) 1 else 0)
chronic.bi02[is.na(clinical.cgvhd)] <- NA
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

pdf('./results/count_histograms.pdf', height=7, width=7)
ggplot(tmp.stats, aes(x=Count)) +
  geom_histogram() +
  facet_wrap(tmp.stats$Measurement %>% factor, scales='free', labeller=as_labeller(tmp.stats.names)) + 
  xlab('number of peptides') +
  theme_light() +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1))
dev.off()


# plot estimates for cGvHD
p2 <- ggplot(tmp %>% na.omit, aes(x=cGvHD %>% factor(levels=c('No','Limited', 'Extensive'), ordered=T), y=Count)) +
  geom_boxplot(outlier.shape=NA, color='grey70') + 
  geom_point(size=0.9, shape=21, fill='grey20', color='grey10', alpha=0.35,  position=position_jitter(0.15)) +
  facet_wrap(~Ligand, scales='free_x', labeller=as_labeller(facet.labs), strip.position='top') +
  theme_light() +
  xlab('cGvHD gradus') +
  ylab('Ligand count') +
  ggtitle('')
  

## Regression models

tr.year <- str_split_fixed(clin$Txdate, '\\.', 3)[, 3] %>% 
  str_split_fixed(., '\\(', 2) %>% .[, 1] %>% 
  str_split_fixed(., '\\.', 5) %>% .[, 1] %>% as.numeric %>% scale

# chronic GvHD
all.cGvHD.glm <- glm(cGvHD ~ Matches+PepNum+HLAfreqsum+Dir+Age+HLAs+Year, data=data.frame(cGvHD=chronic.bi03, 
                                                                                     iedb.matches.allpep %>% dplyr::select(., -ID) %>% scale,
                                                                                     Age=clin$dage %>% scale, 
                                                                                     Dir=tr.dir,
                                                                                     Year=tr.year), family='binomial') %>% summary
imm.cGvHD.glm <- glm(cGvHD ~ Matches+PepNum+HLAfreqsum+Dir+Age+HLAs+Year, data=data.frame(cGvHD=chronic.bi03, 
                                                                                     iedb.matches.immpep %>% dplyr::select(., -ID) %>% scale,
                                                                                     Age=clin$dage %>% scale, 
                                                                                     Dir=tr.dir,
                                                                                     Year=tr.year), family='binomial') %>% summary
mha.cGvHD.glm <- glm(cGvHD ~ MHAs+PepNum+HLAfreqsum+Dir+Age+HLAs+Year, data=data.frame(cGvHD=chronic.bi03, 
                                                                                  mha.matches.allpep %>% dplyr::select(., -ID) %>% scale,
                                                                                  Age=clin$dage %>% scale,
                                                                                  Dir=tr.dir,
                                                                                  Year=tr.year), family='binomial') %>% summary
aff.cGvHD.glm <- glm(cGvHD ~ RANK+PepNum+HLAfreqsum+Dir+Age+HLAs+Year, data=data.frame(cGvHD=chronic.bi03[aff.ranks.immpep.ind], 
                                                                                  iedb.matches.allpep[aff.ranks.immpep.ind, 3:5] %>% scale,
                                                                                  RANK=aff.ranks.immpep$RANK4 %>% scale, 
                                                                                  Age=clin$dage[aff.ranks.immpep.ind] %>% scale, 
                                                                                  Dir=tr.dir[aff.ranks.immpep.ind],
                                                                                  Year=tr.year[aff.ranks.immpep.ind]), family=binomial) %>% summary
aff.cGvHD.glm.6 <- glm(cGvHD ~ RANK+PepNum+HLAfreqsum+Dir+Age+HLAs+Year, data=data.frame(cGvHD=chronic.bi03[aff.ranks.immpep.ind], 
                                                                                  iedb.matches.allpep[aff.ranks.immpep.ind, 3:5] %>% scale,
                                                                                  RANK=aff.ranks.immpep$RANK6 %>% scale, 
                                                                                  Age=clin$dage[aff.ranks.immpep.ind] %>% scale, 
                                                                                  Dir=tr.dir[aff.ranks.immpep.ind],
                                                                                  Year=tr.year[aff.ranks.immpep.ind]), family=binomial) %>% summary
 

## Collect regression results

processCoef <- function(x) {
  x <- x[!(rownames(x)=='(Intercept)'), ]
  x <- x %>% round(., 3) %>% data.frame
  x[, 1] <- x[, 1] 
  x[, 2] <- x[, 2] 
  x <- data.frame(x, c(x[, 1] + x[, 2]), c(x[, 1] - x[, 2]))
  colnames(x) <- c('estimate', 'str.error', 'z.value', 'p.value', 'plus.error', 'minus.error')
  x <- x[1, ] %>% unlist
  rownames(x) <- NULL
  x
}

cgvhd.coef <- data.frame(
  M1=all.cGvHD.glm$coefficients %>% processCoef,
  M2=imm.cGvHD.glm$coefficients %>% processCoef,
  M3=mha.cGvHD.glm$coefficients %>% processCoef,
  M4=aff.cGvHD.glm$coefficients %>% processCoef,
  COEF= c('estimate', 'str.error', 'z.value', 'p.value', 'plus.error', 'minus.error'),
  Condition='cGvHD') %>% melt

cgvhd.coef <- cgvhd.coef %>% spread(., COEF, value) %>% mutate(., BH=p.adjust(p.value, method='BH'))

p4 <- ggplot(cgvhd.coef, aes(y=estimate, x=variable, ymin=minus.error, ymax=plus.error)) +
  geom_pointrange(color='grey1', fatten=8, size=0.25) +
  geom_hline(yintercept=0, color='grey', size=0.4) +
  xlab('') + ylab('Regression coefficient') +
  theme_light() +
  ggtitle('') +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

p5 <- ggplot(cgvhd.coef, aes(y=BH %>% -log10(.), x=variable)) +
  geom_bar(fill='grey50', alpha=0.7, stat='identity') +
  geom_hline(yintercept=-log10(0.05), color='grey', size=0.4) +
  xlab('Analysis method') + ylab(expression('-log'[10]*'(FDR)')) +
  theme_light() +
  theme(strip.background=element_blank(), strip.text.x=element_blank())

pdf(height=4.5, width=7, './results/Figure2.pdf')
ggarrange(p2, 
          ggarrange(p4, p5, nrow=2, ncol=1, align='v', heights=c(1, 1)),
          nrow=1, ncol=2, align='h', labels=c('A', 'B'), widths=c(1, 0.7))
dev.off()

