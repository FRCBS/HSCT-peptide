
## ---------------------------------------------
## Analysis of ligand count data
## ---------------------------------------------

library(data.table)
library(tidyverse)
library(ggpubr)
library(scales)
library(MASS)

setwd("./data/vcf/proteome")



## Ligand count data

iedb.matches.allpep <- fread("./data/IEDB_Ligand_matches.txt", data.table=F)
iedb.matches.allpep <- iedb.matches.allpep[!iedb.matches.allpep$Matches==0, ]
iedb.matches.allpep <- dplyr::rename(iedb.matches.allpep, ID=Sample)

iedb.matches.immpep <- fread("./data/IEDB_Ligand_matches_imm.txt", data.table=F)
iedb.matches.immpep <- dplyr::rename(iedb.matches.immpep, ID=Sample)

mha.matches.allpep  <- fread("./data/IEDB_mHA_matches.txt", data.table=F)
mha.matches.allpep  <- dplyr::rename(mha.matches.allpep, ID=Sample)

aff.ranks.immpep    <- fread("./data/Ranks.tsv", data.table=F)
aff.ranks.immpep$ID <- gsub('T', 'DT', aff.ranks.immpep$ID)


## HLA matching over 6 genes

hlamatch <- fread('~/Projects/HLA_imputation_II/data/McGill_HLAmatch', data.table=F, col.names=c('ID', 'HLAmatch'))


## Parse clinical data

clin    <- read.delim('./data/kaikki_kliiniset_tiedot_1216.txt', stringsAsFactors=F)
clin$ID <- gsub('T', 'DT', clin$ID)
clin    <- clin[match(iedb.matches.allpep$ID, clin$ID), ]
clin    <- clin[!is.na(clin$ID), ]

iedb.matches.allpep <- iedb.matches.allpep[match(clin$ID, iedb.matches.allpep$ID) , ]
iedb.matches.immpep <- iedb.matches.immpep[match(clin$ID, iedb.matches.immpep$ID), ]
mha.matches.allpep  <- mha.matches.allpep[match(clin$ID, mha.matches.allpep$ID), ]
mha.matches.immpep  <- mha.matches.immpep[match(clin$ID, mha.matches.immpep$ID), ]
aff.ranks.immpep    <- aff.ranks.immpep[match(clin$ID, aff.ranks.immpep$ID), ]

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

clin <- left_join(clin, hlamatch, by='ID')


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
chronic.bi03 <- as.numeric(chronic.bi03)


## Summary stats

tmp.stats <- c(AllPeptides_=iedb.matches.allpep$PepNum,
               ImmPeptides_=iedb.matches.immpep$PepNum,
               AllPeptidesIEDB_=iedb.matches.allpep$Matches,
               ImmPeptidesIEDB_=iedb.matches.immpep$Matches,
               AllPeptidesmHA_=mha.matches.allpep$MHAs,
               HighAffinity=aff.ranks.immpep$RANK4)
tmp.stats <- data.frame(Measurement=names(tmp.stats) %>% str_split_fixed(., pattern='_', 3) %>% .[, 1] ,
                        Count=tmp.stats, stringsAsFactors=F)

tmp.stats.names <- c(`AllPeptides`    = "All peptides",
                     `AllPeptidesIEDB`= "All peptides & IEDB",
                     `AllPeptidesmHA` = "All peptides & mHAs",
                     `HighAffinity`   = "Predicted high affinity", 
                     `ImmPeptides`    = "Predicted immunogenic peptides",
                     `ImmPeptidesIEDB`= "Predicted immunogenic peptides\n& IEDB")


## Ligands counts in clinical outcomes

tmp <- full_join(iedb.matches.allpep, clin)[, c('ID', 'Matches', 'cGvHD')] %>% 
  dplyr::rename(., All_ligands=Matches) %>% 
  full_join(., iedb.matches.immpep[, c('ID', 'Matches')]) %>% dplyr::rename(., Immunogenic_ligands=Matches) %>% 
  full_join(., mha.matches.allpep[, c('ID', 'MHAs')]) %>% dplyr::rename(., mHAs=MHAs) %>% 
  full_join(., aff.ranks.immpep[, c('ID', 'RANK4')]) %>% dplyr::rename(., Rank=RANK4) %>% 
  dplyr::select(., -ID)

tmp <- gather(tmp, Ligand, Count, c(1, 3:5))

facet.labs <- c(`All_ligands`="M1", # All HLA ligands
                `Immunogenic_ligands`="M2", # Immunogenic HLA ligands
                `mHAs`="M3", # minor H antigens
                `Rank`="M4") # Predicted affinity rank < 4


# cGvHD ligand count plot
p2 <- ggplot(tmp %>% na.omit, aes(x=cGvHD %>% factor(levels=c('No','Limited', 'Extensive'), ordered=T), y=Count)) +
  geom_boxplot(outlier.shape=NA, color='grey70') + 
  geom_point(size=0.9, shape=1, fill='grey20', color='grey10', alpha=0.9,  position=position_jitter(0.15)) +
  facet_wrap(~Ligand, scales='free', labeller=as_labeller(facet.labs), strip.position='top') +
  xlab('cGvHD gradus') +
  ylab('Ligand count') +
  ggtitle('') +
  theme_minimal() +
  theme(axis.line = element_line(colour="black", size=.2),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(color='black'),
        strip.text=element_text(face='bold')) 


## Regression models

# process transplantation year
tr.year <- str_split_fixed(clin$Txdate, '\\.', 3)[, 3] %>% 
  str_split_fixed(., '\\(', 2) %>% .[, 1] %>% 
  str_split_fixed(., '\\.', 5) %>% .[, 1] %>% as.numeric %>% scale

# chronic GvHD
all.cGvHD.glm <- glm(cGvHD ~ Matches+PepNum+HLAfreqsum+Dir+Age+HLAs+Year+HLAmatch, 
                     data=data.frame(cGvHD=chronic.bi03, 
                                     iedb.matches.allpep %>% dplyr::select(., -ID) %>% scale,
                                     Age=clin$dage %>% scale, 
                                     HLAmatch=clin$HLAmatch,
                                     Dir=tr.dir,
                                     Year=tr.year), family='binomial') 
imm.cGvHD.glm <- glm(cGvHD ~ Matches+PepNum+HLAfreqsum+Dir+Age+HLAs+Year+HLAmatch, 
                     data=data.frame(cGvHD=chronic.bi03, 
                                     iedb.matches.immpep %>% dplyr::select(., -ID) %>% scale,
                                     Age=clin$dage %>% scale,
                                     HLAmatch=clin$HLAmatch,
                                     Dir=tr.dir,
                                     Year=tr.year), family='binomial') 
mha.cGvHD.glm <- glm(cGvHD ~ MHAs+PepNum+HLAfreqsum+Dir+Age+HLAs+Year+HLAmatch, 
                     data=data.frame(cGvHD=chronic.bi03, 
                                     mha.matches.allpep %>% dplyr::select(., -ID) %>% scale,
                                     Age=clin$dage %>% scale,
                                     HLAmatch=clin$HLAmatch,
                                     Dir=tr.dir,
                                     Year=tr.year), family='binomial') 
aff.cGvHD.glm <- glm(cGvHD ~ RANK+PepNum+HLAfreqsum+Dir+Age+HLAs+Year+HLAmatch, 
                     data=data.frame(cGvHD=chronic.bi03, 
                                     iedb.matches.immpep[, 3:5] %>% scale,
                                     RANK=aff.ranks.immpep$RANK4 %>% scale, 
                                     Age=clin$dage %>% scale,
                                     HLAmatch=clin$HLAmatch,
                                     Dir=tr.dir,
                                     Year=tr.year), family=binomial) 
aff.cGvHD.glm.6 <- glm(cGvHD ~ RANK+PepNum+HLAfreqsum+Dir+Age+HLAs+Year+HLAmatch, 
                       data=data.frame(cGvHD=chronic.bi03, 
                                       iedb.matches.allpep[, 3:5] %>% scale,
                                       RANK=aff.ranks.immpep$RANK6 %>% scale, 
                                       Age=clin$dage%>% scale, 
                                       HLAmatch=clin$HLAmatch,
                                       Dir=tr.dir,
                                       Year=tr.year), family=binomial) 

## Collect regression results

# tidy up coefficients table
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
  M1=all.cGvHD.glm %>% summary %>% .$coefficients %>% processCoef,
  M2=imm.cGvHD.glm %>% summary %>% .$coefficients %>% processCoef,
  M3=mha.cGvHD.glm %>% summary %>% .$coefficients %>% processCoef,
  M4=aff.cGvHD.glm %>% summary %>% .$coefficients %>% processCoef,
  COEF= c('estimate', 'str.error', 'z.value', 'p.value', 'plus.error', 'minus.error'),
  Condition='cGvHD') %>% melt

cgvhd.coef <- cgvhd.coef %>% spread(., COEF, value) %>% mutate(., BH=p.adjust(p.value, method='BH'))

# coef + conf. interval plot
p4 <- ggplot(cgvhd.coef, aes(y=estimate, x=variable, ymin=minus.error, ymax=plus.error)) +
  geom_pointrange(color='grey1', fatten=8, size=0.25) +
  geom_hline(yintercept=0, color='grey20', size=0.2, linetype='dashed') +
  xlab('') + ylab('Regression coefficient') +
  ggtitle('') +
  theme_minimal() +
  theme(axis.line = element_line(colour="black", size=.2),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(color='black'),
        strip.text=element_text(face='bold')) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# -log10 p plot
p5 <- ggplot(cgvhd.coef, aes(y=p.value %>% log10 %>% `*`(-1), x=variable)) +
  geom_bar(fill='grey60', alpha=0.7, stat='identity', width=0.5) +
  scale_y_continuous(expand=c(0, 0)) +
  geom_hline(yintercept= 0.038 %>% log10 %>% `*`(-1), color='grey20', size=0.2, linetype='dashed') +
  geom_text(x=3.85, y=1.55, label='FDR <0.05', colour='grey50', size=2.5) +
  xlab('Analysis method') + ylab(expression('-log'[10]*'(p)')) +
  theme_minimal() +
  theme(strip.background=element_blank(), strip.text.x=element_blank()) +
  theme(axis.line = element_line(colour="black", size=.2),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(color='black'),
        strip.text=element_text(face='bold')) 


# collect full regression results

formatCoefTable <- function(x) {
  x <- coefficients(summary(x))
  rownames(x) <- c('(Intercept)', 'Ligand count', 'Number of mismatched peptides', 'Sum of HLA frequencies',
                   'Tr direction', 'Donor age', 'Number of unique HLAs', 'Tr year', 'HLA matching')
  colnames(x) <- c('Estimate', 'Standard deviation', 'z-value', 'p-value')
  x <- signif(x, 3)
  data.frame(Variable=rownames(x), x)
}

# write full regression results
formatCoefTable(all.cGvHD.glm) %>% write.table(., './results/coeff_M1', sep='\t', row.names=F)
formatCoefTable(imm.cGvHD.glm) %>% write.table(., './results/coeff_M2', sep='\t', row.names=F)
formatCoefTable(mha.cGvHD.glm) %>% write.table(., './results/coeff_M3', sep='\t', row.names=F)
formatCoefTable(aff.cGvHD.glm) %>% write.table(., './results/coeff_M4', sep='\t', row.names=F)


# Logistic regression model predictions as probabilities

tmp.dat <- data.frame(iedb.matches.allpep %>% dplyr::select(., -ID) %>% scale,
                      Age=clin$dage %>% scale, 
                      HLAmatch=clin$HLAmatch,
                      Dir=tr.dir,
                      Year=tr.year)
m1.probs <- data.frame(LCount=iedb.matches.allpep$Matches, 
                       Prob=predict(all.cGvHD.glm, newdata=tmp.dat, type='response'),
                       Method='M1')
m1.probs <- data.frame(m1.probs, 
                       predict(lm(Prob~splines::bs(LCount, 3), data=m1.probs), m1.probs, interval='prediction'))

tmp.dat <- data.frame(iedb.matches.immpep %>% dplyr::select(., -ID) %>% scale,
                      Age=clin$dage %>% scale,
                      HLAmatch=clin$HLAmatch,
                      Dir=tr.dir,
                      Year=tr.year)
m2.probs <- data.frame(LCount=iedb.matches.immpep$Matches, 
                       Prob=predict(imm.cGvHD.glm, newdata=tmp.dat, type='response'),
                       Method='M2')
m2.probs <- data.frame(m2.probs, 
                       predict(lm(Prob~splines::bs(LCount, 3), data=m2.probs), m2.probs, interval='prediction'))

tmp.dat <- data.frame(mha.matches.allpep %>% dplyr::select(., -ID) %>% scale,
                      Age=clin$dage %>% scale,
                      HLAmatch=clin$HLAmatch,
                      Dir=tr.dir,
                      Year=tr.year)
m3.probs <- data.frame(LCount=mha.matches.allpep$MHAs, 
                       Prob=predict(mha.cGvHD.glm, newdata=tmp.dat, type='response'),
                       Method='M3')
m3.probs <- data.frame(m3.probs, 
                       predict(lm(Prob~splines::bs(LCount, 3), data=m3.probs), m3.probs, interval='prediction'))

tmp.dat <- data.frame(iedb.matches.allpep[, 3:5] %>% scale,
                      RANK=aff.ranks.immpep$RANK4 %>% scale, 
                      Age=clin$dage %>% scale,
                      HLAmatch=clin$HLAmatch,
                      Dir=tr.dir,
                      Year=tr.year)
m4.probs <- data.frame(LCount=aff.ranks.immpep$RANK4, 
                       Prob=predict(aff.cGvHD.glm, newdata=tmp.dat, type='response'),
                       Method='M4')
m4.probs <- data.frame(m4.probs, 
                       predict(lm(Prob~splines::bs(LCount, 3), data=m4.probs), m4.probs, interval='prediction'))

# severe cGvHD probability vs. ligand count plots
p6 <- ggplot(rbind(m1.probs, m2.probs, m3.probs, m4.probs), aes(LCount, Prob, ymin=lwr, ymax=upr)) +
  geom_point(size=0.9, shape=1, fill='grey20', color='grey10', alpha=0.9) +
  geom_line(aes(LCount, fit), rbind(m1.probs, m2.probs, m3.probs, m4.probs), size=0.2) +
  geom_ribbon(alpha=0.5, fill='grey85') +
  ggtitle('') +
  theme_minimal() +
  xlab('Ligand count') + ylab('Probability of severe cGvHD') +
  facet_wrap(~Method, scales='free_x', ncol=4) +
  theme(strip.background=element_blank()) +
  theme(axis.line = element_line(colour="black", size=0),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(color='black'),
        strip.text=element_text(face='bold')) +
  coord_cartesian(ylim=c(0, 1)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=.2) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=.2)



## Collect all plots  

pdf(height=6.8, width=7, './results/Figure3_Ligand_effects.pdf')
ggarrange(ggarrange(p2, 
                    ggplot() + theme_void(),
                    ggarrange(p4, p5, nrow=2, ncol=1, align='v', heights=c(1, 1)),
                    nrow=1, ncol=3, align='h', labels=c('A', '', 'B'), widths=c(1, 0.05, 0.55)),
          ggarrange(p6, labels='C'), 
          nrow=2, ncol=1, heights=c(1, 0.55))
dev.off()


