library(readr)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(Rmisc)
library(ez)
library(psychReport)
library(singcar)
library(tidyverse)
library(Hmisc)

#PATHS: will need changing
# path for data & analysis
anaPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/analysis/radial_reaching'

# load data
setwd(anaPath)
res <- read.csv('all_radial-reaching_compiled.csv')

# summary data
# extracting data from furthest two target locs
res_periph <- subset(res, res$POSITION == -400 | res$POSITION == -300 | 
                         res$POSITION == 300 | res$POSITION == 400)
res_medians <- aggregate(AE ~ PPT*POSITION*VIEW*SIDE*DOM*DIAGNOSIS*GRP*SITE*AGE, 
                         median, data = res_periph)
colnames(res_medians)[colnames(res_medians)=='AE'] <- 'AEmed'

# removing free + peripheral trails at 100mm left for 101, error
res_medians <- res_medians[order(res_medians$PPT) ,]

res_means <- aggregate(AEmed ~ PPT*VIEW*SIDE*DOM*GRP*DIAGNOSIS*SITE*AGE, 
                       mean, data = res_medians)
colnames(res_means)[colnames(res_means)=='AEmed'] <- 'AEmean'

# casting by task
PMIdata <- dcast(res_means, PPT+GRP+DIAGNOSIS+DOM+SIDE+SITE+AGE ~ VIEW)
PMIdata$PMI <- PMIdata$Peripheral - PMIdata$Free
PMIgrand <- aggregate(PMI ~ PPT*GRP*DIAGNOSIS*SITE*AGE, mean, data = PMIdata)

write.csv(res_medians, 'radial-medians.csv', row.names = FALSE)
write.csv(res_means, 'radial-means.csv', row.names = FALSE)
write.csv(PMIdata, 'radialPMI.csv', row.names = FALSE)

## summary data
meanPMI_side <- summarySE(PMIdata, measurevar = 'PMI', groupvar = c('DIAGNOSIS', 'DOM'),
                          na.rm = TRUE)
write.csv(meanPMI_side, 'radialPMI_means.csv', row.names = FALSE)

## getting & saving demographics for pub
AGE <- summarySE(PMIgrand, measurevar = 'AGE', groupvar = c('DIAGNOSIS'), na.rm = TRUE)
res$ED <- as.numeric(res$ED)
ED <- aggregate(ED ~ PPT * GRP * DIAGNOSIS, mean, data = res)
ED <- summarySE(ED, measurevar = 'ED', groupvar = c('DIAGNOSIS'), na.rm = TRUE)

##### ANOVA - differences between sites #####
siteANOVA <- na.omit(res_means)
siteANOVA <- siteANOVA[siteANOVA$PPT != 212 ,]
siteANOVA <- siteANOVA[siteANOVA$PPT != 310 , c(1:7,9)]
siteANOVA <- aggregate(AEmean ~ PPT+VIEW+DOM+DIAGNOSIS+SITE, mean, data = siteANOVA)

SITE_ANOVA <- ezANOVA(
  data = siteANOVA
  , dv = .(AEmean)
  , wid = .(PPT)
  , within = .(VIEW)
  , between =. (SITE)
  , type = 3
)

print(SITE_ANOVA)
# no difference across sites! Hooray!

#### outlier removal step has been removed from this analysis - code with this in can be found:
#### DMT -> analysis -> radial_reaching -> as-preregistered

######## CASE CONTROL ANALYSIS ########
# create data-frames (controls and patients) with key information
td_summary <- summarySE(data = PMIdata, measurevar = 'PMI', 
                        groupvars = c('DIAGNOSIS','DOM','SITE'), na.rm = TRUE)
age <- summarySE(data = PMIdata, measurevar = 'AGE', 
                 groupvars = c('DIAGNOSIS','DOM','SITE'), na.rm = TRUE)
age <- age[, c(1:3,5)]
td_summary <- merge(td_summary, age)

#patient data for test of deficit
td_patient <- PMIdata[, c(1,3,4,6,7,10)]
td_patient <- td_patient[td_patient$DIAGNOSIS != 'HC' ,] #removing controls
# next, run case-control independently for each site - different set ups

## UOE CASE CONTROL ##
## correlation matrix of age and PMI - for BTD_cov
## non-dominant
PMIcorr <- na.omit(PMIdata)
PMIcorr <- PMIcorr[PMIcorr$DOM == 'ND' & PMIcorr$SITE == 'UOE', c(10,7)]
PMIcorr$AGE <- as.numeric(PMIcorr$AGE)
corr <- rcorr(as.matrix(PMIcorr))
# relabelling, non-dominant correlation matrix
NDCM_UOE <- corr[["r"]]

## dominant
PMIcorr <- na.omit(PMIdata)
PMIcorr <- PMIcorr[PMIcorr$DOM == 'D' & PMIcorr$SITE == 'UOE', c(10,7)]
PMIcorr$AGE <- as.numeric(PMIcorr$AGE)
corr <- rcorr(as.matrix(PMIcorr))
# relabelling, non-dominant correlation matrix
DCM_UOE <- corr[["r"]]

# isolating control data for analysis
tdUOE_control <- td_summary[td_summary$DIAGNOSIS == 'HC' 
                            & td_summary$SITE == 'UOE' ,]
# matrices for input into BTD_cov function
tdUOE_controlND <- tdUOE_control[tdUOE_control$DOM == 'ND', c(5,6)]
tdUOE_controlD <- tdUOE_control[tdUOE_control$DOM == 'D', c(5,6)]

# patient data
tdUOE_patient <- td_patient[td_patient$SITE == 'UOE' ,]
tdUOE_patient <- dcast(tdUOE_patient, PPT+DIAGNOSIS+AGE ~ DOM)
# NA values  = -1, so t-value is negative and can remove later
tdUOE_patient[is.na(tdUOE_patient$D), "D"] <- -1 

# time for test of deficit! Calling on 'singcar' package developed by Jonathan Rittmo
# using Crawford et al.'s (2011) Bayesian test of covariates (with singcar package)
td_dom <- read.csv(text = 'PMI,TSTAT,PVALUE,PROP,ZCC,CI,LCI-T,HCI-T,LCI-P,HCI-P,DOM,PPT,DIAGNOSIS')
td_ndom <- read.csv(text = 'PMI,TSTAT,PVALUE,PROP,ZCC,CI,LCI-T,HCI-T,LCI-P,HCI-P,DOM,PPT,DIAGNOSIS')
for (l in 1:length(tdUOE_patient$PPT)){
  #left data first
  NDres <- BTD_cov(tdUOE_patient$ND[l], tdUOE_patient$AGE[l], tdUOE_controlND, tdUOE_control$AGE[1], 
                   alternative = 'greater', int_level = 0.95, iter = 10000,
                   use_sumstats = TRUE, cor_mat = NDCM_UOE, sample_size = tdUOE_control$N[1])
  diff <- t(NDres$estimate)
  ltmp <- data.frame(tdUOE_patient$ND[l], NDres$statistic, NDres$p.value, 
                     diff[1], diff[2], t(NDres$interval), 'ND', check.names = FALSE) 
  ltmp$PPT <- tdUOE_patient$PPT[l]
  ltmp$DIAGNOSIS <- tdUOE_patient$DIAGNOSIS[l]
  td_ndom <- rbind(td_ndom, ltmp)
  #then right data
  Dres <- BTD_cov(tdUOE_patient$D[l], tdUOE_patient$AGE[l], tdUOE_controlD, tdUOE_control$AGE[2], 
                  alternative = 'greater', int_level = 0.95, iter = 10000,
                  use_sumstats = TRUE, cor_mat = DCM_UOE, sample_size = tdUOE_control$N[2])
  diff <- t(Dres$estimate)
  rtmp <- data.frame(tdUOE_patient$D[l], Dres$statistic, Dres$p.value, 
                     diff[1], diff[2], t(Dres$interval), 'D', check.names = FALSE) 
  rtmp$PPT <- tdUOE_patient$PPT[l]
  rtmp$DIAGNOSIS <- tdUOE_patient$DIAGNOSIS[l]
  td_dom <- rbind(td_dom, rtmp)
}

# merging and renaming data-frames
tdUOE_results <- read.csv(text = 'PMI,TSTAT,PVALUE,ZCC,PROP,CI,LCI-T,HCI-T,LCI-P,HCI-P,DOM,PPT,DIAGNOSIS')
# changing names of td data-frames to match td-res
names(td_ndom) <- names(tdUOE_results)
names(td_dom) <- names(tdUOE_results)
tdUOE_results <- rbind(tdUOE_results, td_ndom, td_dom)
# remove original NA values (PMI = -1), where patients had no data
tdUOE_results <- tdUOE_results[tdUOE_results$PMI > -0.9 ,]
tdUOE_results$PPT <- factor(tdUOE_results$PPT)
# add site back to results
tdUOE_results$SITE <- 'UOE'

# identifying patients with significant deficit
tdUOE_results$DEFICIT <- tdUOE_results$PVALUE < 0.025
# identifying patients with possible deficit
tdUOE_results$BL <- tdUOE_results$PVALUE < 0.05
# convert logicals into numerics, useful for binomial later
tdUOE_results$DEFICIT <- as.numeric(tdUOE_results$DEFICIT)
tdUOE_results$BL <- as.numeric(tdUOE_results$BL)

## UEA CASE CONTROL ##
## correlation matrix of age and PMI - for BTD_cov
## non-dominant
PMIcorr <- na.omit(PMIdata)
PMIcorr <- PMIcorr[PMIcorr$DOM == 'ND' & PMIcorr$SITE == 'UEA', c(10,7)]
PMIcorr$AGE <- as.numeric(PMIcorr$AGE)
corr <- rcorr(as.matrix(PMIcorr))
# relabelling, non-dominant correlation matrix
NDCM_UEA <- corr[["r"]]

## dominant
PMIcorr <- na.omit(PMIdata)
PMIcorr <- PMIcorr[PMIcorr$DOM == 'D' & PMIcorr$SITE == 'UEA', c(10,7)]
PMIcorr$AGE <- as.numeric(PMIcorr$AGE)
corr <- rcorr(as.matrix(PMIcorr))
# relabelling, non-dominant correlation matrix
DCM_UEA <- corr[["r"]]

# isolating control data for analysis
tdUEA_control <- td_summary[td_summary$DIAGNOSIS == 'HC' 
                            & td_summary$SITE == 'UEA' ,]
# matrices for input into BTD_cov function
tdUEA_controlND <- tdUEA_control[tdUEA_control$DOM == 'ND', c(5,6)]
tdUEA_controlD <- tdUEA_control[tdUEA_control$DOM == 'D', c(5,6)]

#patient data
tdUEA_patient <- td_patient[td_patient$SITE == 'UEA' ,]
tdUEA_patient <- dcast(tdUEA_patient, PPT+DIAGNOSIS+AGE ~ DOM)

# time for test of deficit! Calling on 'singcar' package developed by Jonathan Rittmo
# using Crawford et al's (2011) Bayesian test of covariates (with singcar package)
td_dom <- read.csv(text = 'PMI,TSTAT,PVALUE,PROP,ZCC,CI,LCI-T,HCI-T,LCI-P,HCI-P,DOM,PPT,DIAGNOSIS')
td_ndom <- read.csv(text = 'PMI,TSTAT,PVALUE,PROP,ZCC,CI,LCI-T,HCI-T,LCI-P,HCI-P,DOM,PPT,DIAGNOSIS')
for (l in 1:length(tdUEA_patient$PPT)){
  #left data first
  NDres <- BTD_cov(tdUEA_patient$ND[l], tdUEA_patient$AGE[l], tdUEA_controlND, tdUEA_control$AGE[1], 
              alternative = 'greater', int_level = 0.95, iter = 10000,
              use_sumstats = TRUE, cor_mat = NDCM_UEA, sample_size = tdUEA_control$N[1])
  diff <- t(NDres$estimate)
  ltmp <- data.frame(tdUEA_patient$ND[l], NDres$statistic, NDres$p.value, 
                     diff[1], diff[2], t(NDres$interval), 'ND', check.names = FALSE) 
  ltmp$PPT <- tdUEA_patient$PPT[l]
  ltmp$DIAGNOSIS <- tdUEA_patient$DIAGNOSIS[l]
  td_ndom <- rbind(td_ndom, ltmp)
  #then right data
  Dres <- BTD_cov(tdUEA_patient$D[l], tdUEA_patient$AGE[l], tdUEA_controlD, tdUEA_control$AGE[2], 
                  alternative = 'greater', int_level = 0.95, iter = 10000,
                  use_sumstats = TRUE, cor_mat = DCM_UEA, sample_size = tdUEA_control$N[2])
  diff <- t(Dres$estimate)
  rtmp <- data.frame(tdUEA_patient$D[l], Dres$statistic, Dres$p.value, 
                     diff[1], diff[2], t(Dres$interval), 'D', check.names = FALSE) 
  rtmp$PPT <- tdUEA_patient$PPT[l]
  rtmp$DIAGNOSIS <- tdUEA_patient$DIAGNOSIS[l]
  td_dom <- rbind(td_dom, rtmp)
}

# merging and renaming data-frames
tdUEA_results <- read.csv(text = 'PMI,TSTAT,PVALUE,ZCC,PROP,CI,LCI-T,HCI-T,LCI-P,HCI-P,DOM,PPT,DIAGNOSIS')
# changing names of td data-frames to match td-res
names(td_ndom) <- names(tdUEA_results)
names(td_dom) <- names(tdUEA_results)
tdUEA_results <- rbind(tdUEA_results, td_ndom, td_dom)
# remove original NA values (PMI = -1), where patients had no data
tdUEA_results <- tdUEA_results[tdUEA_results$PMI > -0.9 ,]
tdUEA_results$PPT <- factor(tdUEA_results$PPT)
# add site back to results
tdUEA_results$SITE <- 'UEA'

# identifying patients with significant deficit
tdUEA_results$DEFICIT <- tdUEA_results$PVALUE < 0.025
# identifying patients with possible deficit
tdUEA_results$BL <- tdUEA_results$PVALUE < 0.05
# convert logicals into numerics, useful for binomial later
tdUEA_results$DEFICIT <- as.numeric(tdUEA_results$DEFICIT)
tdUEA_results$BL <- as.numeric(tdUEA_results$BL)

## merge UOE and UEA data-frames
td_results <- rbind(tdUOE_results, tdUEA_results)
#save
write.csv(td_results, 'radial-reaching_case-control.csv', row.names = FALSE)

##### BINOMIAL STATS #####
# calculating the likelihood that inflated PMI occurs at above chance in each patient group
# with filtered data. Deficit and borderline deficit
pval <- .05
pval_Bl <- .1

## take any p < .025 for either side, re-run binomial on this data
# extracting p-value for each side (left and right) using td_results
td_side <- dcast(td_results, PPT+DIAGNOSIS+SITE ~ DOM, value.var = 'PVALUE')
td_side[is.na(td_side$D), "D"] <- 10
# finding defcits and borderline cases
td_side$DEFICIT <- (td_side$ND < 0.025 | td_side$D < 0.025)
td_side$BL <- (td_side$ND < 0.05 | td_side$D < 0.05)
# changing to numerical value
td_side$DEFICIT <- as.numeric(td_side$DEFICIT)
td_side$BL <- as.numeric(td_side$BL)

# splitting into MCI and AD
MCIP <- td_side[td_side$DIAGNOSIS == 'MCI', ]
ADP <- td_side[td_side$DIAGNOSIS == 'AD', ]

# binomial test of deficit
binMCI <- binom.test(sum(MCIP$DEFICIT), length(MCIP$PPT), pval, alternative = 'greater')
print(binMCI)
binAD <- binom.test(sum(ADP$DEFICIT), length(ADP$PPT), pval, alternative = 'greater')
print(binAD)
binALL <- binom.test(sum(td_side$DEFICIT), length(td_side$PPT), pval, alternative = 'greater')
print(binALL)

# binomial test of borderline deficit
binMCI <- binom.test(sum(MCIP$BL), length(MCIP$PPT), pval_Bl, alternative = 'greater')
print(binMCI)
binAD <- binom.test(sum(ADP$BL), length(ADP$PPT), pval_Bl, alternative = 'greater')
print(binAD)
binALL <- binom.test(sum(td_side$BL), length(td_side$PPT), pval_Bl, alternative = 'greater')
print(binALL)

##### ANOVAS, group comparisons #####
# use PMI data-frame to run between-subject ANOVA
# collapse across side to run anova
PMIanova <- aggregate(PMI ~ PPT+DIAGNOSIS+SITE+AGE, mean, data = PMIdata)

# FULL ANOVA ON FILTERED DATA
PMI_ANOVA <- ezANOVA(
  data = PMIanova
  , dv = .(PMI)
  , wid = .(PPT)
  , between = .(DIAGNOSIS)
  , between_covariates = .(SITE, AGE)
  , type = 3,
  return_aov = TRUE,
  detailed = TRUE
)

PMI_ANOVA$ANOVA
PMI_ANOVA$`Mauchly's Test for Sphericity`
PMI_ANOVA$`Sphericity Corrections`
aovPMI <- aovEffectSize(ezObj = PMI_ANOVA, effectSize = "pes")
aovDispTable(aovPMI)

#pair-wise t-test
PMIttest <- pairwise.t.test(PMIanova$PMI, PMIanova$DIAGNOSIS, p.adj = 'bonf')
print(PMIttest)

##### PLOTTING #####
## PLOT: PMI side + av ##
# make control data-frame
control_PMI <- subset(PMIdata, PMIdata$DIAGNOSIS == 'HC')
control_PMI$TSTAT <- 0
control_PMI$PVALUE <- 1
control_PMI$DEFICIT <- 0
control_PMI$BL <- 0
# get deficit data for patients
plot_PMI <- merge(PMIdata, td_results, by = c('PPT','DOM','DIAGNOSIS','SITE'))
# include only relevant info
plot_PMI$PMI <- plot_PMI$PMI.x
plot_PMI <- plot_PMI[, c(1:9,23,12,13,21,22)]

# make plot data frame
plot_PMI <- rbind(control_PMI, plot_PMI)

plot_PMI$DOM <- factor(plot_PMI$DOM, levels = c('ND', 'D'))
plot_PMI$DIAGNOSIS <- factor(plot_PMI$DIAGNOSIS, levels = c('HC','MCI','AD'))
plot_PMI$PPT <- factor(plot_PMI$PPT)
# make column where deficit and BL cases are combined
# 1 = no deficit, 2 = borderline deficit, 3 = deficit
plot_PMI$DEFICITS <- as.numeric(plot_PMI$DEFICIT) + as.numeric(plot_PMI$BL)
plot_PMI$DEFICITS <- factor(plot_PMI$DEFICITS)

ggplot(plot_PMI, aes(x = DOM, y = PMI, colour = DIAGNOSIS, group = PPT, shape = DEFICITS)) + 
  geom_line(aes(group = PPT), alpha = .7, size = 0.7, position = position_dodge(.2)) +
  geom_point(size = 2.5, position = position_dodge(.2)) +
  stat_summary(aes(y = PMI, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1, size = 4, group = 1) +
  scale_color_manual(values = c('grey45','grey45','grey45')) +
  scale_shape_manual(values = c(1,18,18)) +
  facet_wrap(~DIAGNOSIS) +
  labs(x = 'Side', y = 'Radial PMI (mm)') +
  theme_classic() + theme(legend.position = 'none', 
                          axis.title = element_text(size = 12),
                          axis.text = element_text(size = 10),
                          strip.text = element_text(size = 10) 
  ) -> pPMI
pPMI

## PLOT 4: average PMI across sides - combine with PLOT 2
PMIav_plot <- aggregate(PMI ~ PPT*DIAGNOSIS*AGE*SITE, mean, data = plot_PMI)
PMIav_plot <- PMIav_plot[order(PMIav_plot$PPT),]
plot_PMI$DEFICITS <- as.numeric(plot_PMI$DEFICITS)
deficit <- aggregate(DEFICITS ~ PPT*DIAGNOSIS, max, data = plot_PMI)
deficit <- deficit[order(deficit$PPT),]
deficit <- deficit[deficit$PPT != 310 ,]
PMIav_plot$DEFICITS <- factor(deficit$DEFICITS)

ggplot(PMIav_plot, aes(x = DIAGNOSIS, y = PMI, colour = DIAGNOSIS, group = PPT, shape = DEFICITS)) + 
  geom_point(size = 2.5, position = position_dodge(.2)) +
  stat_summary(aes(y = PMI, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1, size = 4, group = 1) +
  scale_shape_manual(values = c(1,18,18)) +
  scale_colour_manual(values = c('grey45','grey45','grey45')) +
  labs(x = 'Diagnosis', y = '') +
  theme_classic() + theme(legend.position = 'none', 
                          axis.title = element_text(size = 12),
                          axis.text = element_text(size = 10),
                          strip.text = element_text(size = 10) 
  ) -> avPMI
avPMI

# combining in to main plot for publication
PMIfig <- ggarrange(pPMI, avPMI,
                    ncol=2, nrow=1,
                    widths = c(1.5,1),
                    labels = c('A','B'),
                    hjust = -1)
PMIfig

ggsave('RADAE-fig.png', plot = last_plot(), device = NULL, dpi = 300, 
       width = 8, height = 4, path = anaPath)
