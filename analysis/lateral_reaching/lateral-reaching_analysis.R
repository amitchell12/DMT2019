##### LATERAL REACHING ANALYSIS CODE
## Code for main analysis for Mitchell et al. ....

library(readr)
library(ggplot2)
library(Rmisc)
library(gridExtra)
library(plyr)
library(tidyverse)
library(reshape2)
library(Hmisc)
library(ggpubr)
library(singcar)
library(ez)
library(psychReport)

#set working directory to where data is
# on mac (desktop)
dataPath <- '/Users/Alex/Documents/DMT/data/'
anaPath <- '/Users/Alex/Documents/DMT/analysis//lateral_reaching/'
# on mac (laptop)
#dataPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/data/'
#anaPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/analysis/lateral_reaching/'
setwd(anaPath)

# load data file
res <- read.csv('lateral-reaching_compiled.csv')

## getting dominant + non-dominant sides, for analysis
names(res)[4] <- 'HAND'
res$HAND <- factor(res$HAND, labels = c('left','right'))
# if hand = side, dominant; else non dominant
res$DOM <- as.numeric(res$SIDE == res$HAND) #1 = dominant, 0 = non-dominant
res$DOM <- factor(res$DOM, labels= c('ND','D'))
# change order so dominance up-front
res <- res[, c(1:8,27,9:26)]

## remove trial with widly high RT (> 6000ms)
res <- res[res$RT < 6000 ,]

# counting total number of valid trials for each patient
valid <- aggregate(TARGx ~ DIAGNOSIS * PPT * DOM * VIEW, length, data = res)
# calculate percentage of trials completed in each group
valid$HITS <- (valid$TARGx/27)*100
valid_tot <- aggregate(HITS ~ DIAGNOSIS * VIEW, mean, data = valid)

######### step 1, calculating PMI, all data ############
res_medians <- aggregate(
  AE ~ PPT * DOM * SIDE * VIEW * POSITION * SITE * GRP * DIAGNOSIS * AGE * ED, 
  median, data = res)
colnames(res_medians)[colnames(res_medians)=='AE'] <- 'AEmed' #change name to be more logical

# changing levels to be more informative
res_medians$GRP <- factor(res_medians$GRP)
levels(res_medians$GRP) <- c('Control', 'Patient') #changing group name from 1 = control, 2 = AD
levels(res_medians$VIEW) <- c('Free', 'Peripheral')
res_medians$SITE <- factor(res_medians$SITE)
levels(res_medians$SITE) <- c('UOE', 'UEA')
res_medians$DIAGNOSIS <- factor(res_medians$DIAGNOSIS)
res_medians <- res_medians[order(res_medians$PPT), ] 

res_means <- aggregate(AEmed ~ PPT * DOM * SIDE * VIEW * SITE * GRP * DIAGNOSIS * AGE * ED, 
                       mean, data = res_medians)
res_means <- res_means[order(res_means$PPT), ] 
colnames(res_means)[colnames(res_means) == 'AEmed'] <- 'AEmean'
# save data
write.csv(res_medians, 'lateral-medians_all.csv', row.names = FALSE)
write.csv(res_means, 'lateral-means_all.csv', row.names = FALSE)

# to calculate PMI need to cast by task....
PMIdata <- dcast(res_means, PPT+GRP+SITE+DOM+SIDE+DIAGNOSIS+AGE+ED ~ VIEW) #different data-frame
PMIdata$PMI <- PMIdata$Peripheral - PMIdata$Free
write.csv(PMIdata, 'lateralPMI_all.csv', row.names = FALSE)

# summary
meanPMI_side <- summarySE(PMIdata, measurevar = 'PMI', groupvar = c('DIAGNOSIS', 'DOM'),
                       na.rm = TRUE)
write.csv(meanPMI_side, 'lateralPMI_means.csv', row.names = FALSE)
# averaged across side
PMIgrand <- aggregate(PMI~PPT+GRP+SITE+DIAGNOSIS+AGE+ED, mean, data=PMIdata)
meanPMI_grand <- summarySE(PMIgrand, measurevar = 'PMI', groupvar = c('DIAGNOSIS'))

## bit of plotting by eccentricity
# controls
res_medians$POSITION <- factor(res_medians$POSITION)
meds_control <- res_medians[res_medians$GRP == 'Control' ,]
  
ggplot(meds_control, aes(x = POSITION, y = AEmed, colour = DOM)) +
  geom_point(shape = 1, size = 2) +
  geom_line(aes(group = PPT), size = 0.5, alpha = .5) +
  facet_grid(cols = vars(DOM), rows = vars(VIEW)) + 
  labs(title = 'Control', x = 'Eccentricity (deg)', 
       y = 'Mean AE (mm)', element_text(size = 12)) +
  scale_colour_manual(values = c('black', 'grey50')) +
  theme_bw() + theme(legend.position = 'none', text = element_text(size = 10),
                    strip.text.x = element_text(size = 10)) 


ggsave('control_ecc.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

# patients - MCI
meds_MCI <- res_medians[res_medians$DIAGNOSIS == 'MCI' ,]

ggplot(meds_MCI, aes(x = POSITION, y = AEmed, colour = DOM)) +
  geom_point(shape = 1, size = 2) +
  geom_line(aes(group = PPT), size = 0.5, alpha = .5) +
  facet_grid(cols = vars(DOM), rows = vars(VIEW)) +
  labs(title = 'MCI', x = 'Eccentricity (deg)', 
       y = 'Mean AE (mm)', element_text(size = 12)) +
  scale_colour_manual(values = c('black', 'grey50')) +
  theme_bw() + theme(legend.position = 'none', text = element_text(size = 10),
                     strip.text.x = element_text(size = 10)) 

ggsave('MCI_ecc.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

# patients - AD
meds_AD <- res_medians[res_medians$DIAGNOSIS == 'AD' ,]

ggplot(meds_AD, aes(x = POSITION, y = AEmed, colour = DOM)) +
  geom_point(shape = 1, size = 2) +
  geom_line(aes(group = PPT), size = 0.5, alpha = .5) +
  facet_grid(cols = vars(DOM), rows = vars(VIEW)) + 
  labs(title = 'Alzheimers', x = 'Eccentricity (deg)', 
       y = 'Mean AE (mm)', element_text(size = 12)) +
  scale_colour_manual(values = c('black', 'grey50')) +
  theme_bw() + theme(legend.position = 'none', text = element_text(size = 10),
                     strip.text.x = element_text(size = 10)) 

ggsave('AD_ecc.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

#### NOTE: outlier removal stage removed from this analysis - can be found in:
#### DMT -> analysis -> lateral_reaching -> as-preregistered

######## step 2: single case stats, filtered data #########
# create data-frames (controls and patients) with key information
td_summary <- summarySE(data = PMIdata, measurevar = 'PMI', 
                        groupvars = c('DIAGNOSIS','DOM'), na.rm = TRUE)
# isolating control data for analysis
td_control <- td_summary[td_summary$DIAGNOSIS == 'HC', ]
#patient data for test of deficit
td_patient <- PMIdata[, c(1,4,6,11)]
td_patient <- td_patient[td_patient$DIAGNOSIS != 'HC' ,] #removing controls
td_patient <- dcast(td_patient, PPT+DIAGNOSIS ~ DOM)
# NA values  = -1, so t-value is negative and can remove later
td_patient[is.na(td_patient$D), "D"] <- -1 

#patient data for test of deficit = same as above
# time for test of deficit! Calling on 'singcar' package developed by Jonathan Rittmo
# using Crawford's (1998) test of deficit
td_dom <- read.csv(text = 'PMI,TSTAT,PVALUE,PROP,ZCC,CI,LCI-T,HCI-T,LCI-P,HCI-P,DOM,PPT,DIAGNOSIS')
td_ndom <- read.csv(text = 'PMI,TSTAT,PVALUE,PROP,ZCC,CI,LCI-T,HCI-T,LCI-P,HCI-P,DOM,PPT,DIAGNOSIS')
for (l in 1:length(td_patient$PPT)){
  #left data first
  NDres <- TD(td_patient$ND[l], td_control$PMI[1], td_control$sd[1], td_control$N[1], 
                alternative = 'greater', na.rm = FALSE)
  diff <- t(NDres$estimate)
  ltmp <- data.frame(td_patient$ND[l], NDres$statistic, NDres$p.value, 
                     diff[1], diff[2], t(NDres$interval), 'ND', check.names = FALSE) 
  ltmp$PPT <- td_patient$PPT[l]
  ltmp$DIAGNOSIS <- td_patient$DIAGNOSIS[l]
  td_ndom <- rbind(td_ndom, ltmp)
  #then right data
  Dres <- TD(td_patient$D[l], td_control$PMI[2], td_control$sd[2], td_control$N[2], 
                 alternative = 'greater', na.rm = FALSE)
  diff <- t(Dres$estimate)
  rtmp <- data.frame(td_patient$D[l], Dres$statistic, Dres$p.value, 
                     diff[1], diff[2], t(Dres$interval), 'D', check.names = FALSE) 
  rtmp$PPT <- td_patient$PPT[l]
  rtmp$DIAGNOSIS <- td_patient$DIAGNOSIS[l]
  td_dom <- rbind(td_dom, rtmp)
}

# merging and renaming data-frames
td_results <- read.csv(text = 'PMI,TSTAT,PVALUE,ZCC,PROP,CI,LCI-T,HCI-T,LCI-P,HCI-P,DOM,PPT,DIAGNOSIS')
# changing names of td data-frames to match td-res
names(td_ndom) <- names(td_results)
names(td_dom) <- names(td_results)
td_results <- rbind(td_results, td_ndom, td_dom)
# remove original NA values (PMI = -1), where patients had no data
td_results <- td_results[td_results$PMI > 0 ,]
td_results$PPT <- factor(td_results$PPT)

# identifying patients with significant deficit
td_results$DEFICIT <- td_results$PVALUE < 0.025
# identifying patients with possible deficit
td_results$BL <- td_results$PVALUE < 0.05
# convert logicals into numerics, useful for binomial later
td_results$DEFICIT <- as.numeric(td_results$DEFICIT)
td_results$BL <- as.numeric(td_results$BL)

#save
write.csv(td_results, 'lateral_case-control.csv', row.names = FALSE)

## plotting p-values
ggplot(td_results, aes(x = DOM, y = PVALUE, group = PPT, colour = DIAGNOSIS)) +
  geom_point(size = 2, position = position_dodge(.2)) +
  geom_line(aes(group = PPT), size = 0.5, alpha = .5, position = position_dodge(.2)) + 
  scale_color_manual(values = c('black', 'grey50')) +
  geom_hline(yintercept = 0.025, color = 'red') +
  geom_hline(yintercept = 0.05, color = 'blue') +
  labs(title = 'Single case stats, filtered', x = 'Side', y = 'Probability of deficit', 
       element_text(size = 10)) +
  theme_bw()

ggsave('SCS_probability-scatter.png', plot = last_plot(), device = NULL, dpi = 300, 
       width = 6, height = 7, path = anaPath)

# Binomial statistics
# calculating the likelihood that inflated PMI occurs at above chance in each patient group
# with filtered data. Deficit and borderline deficit
## take any p < .025 for either side, re-run binomial on this data
# extracting p-value for each side (left and right) using td_results
pval <- .05
pval_Bl <- .1

td_side <- dcast(td_results, PPT+DIAGNOSIS ~ DOM, value.var = 'PVALUE')
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
binAD <- binom.test(sum(ADP$DEFICIT), length(ADP$DEFICIT), pval, alternative = 'greater')
print(binAD)
binALL <- binom.test(sum(td_side$DEFICIT), length(td_side$DEFICIT), pval, alternative = 'greater')
print(binALL)

# binomial test for borderline deficit
binMCI <- binom.test(sum(MCIP$BL), length(MCIP$PPT), pval_Bl, alternative = 'greater')
print(binMCI)
binAD <- binom.test(sum(ADP$BL), length(ADP$PPT), pval_Bl, alternative = 'greater')
print(binAD)
binALL <- binom.test(sum(td_side$BL), length(td_side$PPT), pval_Bl, alternative = 'greater')
print(binALL)

###### ANOVA ######
# FULL ANOVA ON FILTERED DATA
PMI_ANOVA <- ezANOVA(
  data = PMIgrand
  , dv = .(PMI)
  , wid = .(PPT)
  , between = .(DIAGNOSIS)
  , between_covariates = .(AGE)
  , type = 3,
  return_aov = TRUE,
  detailed = TRUE
)

PMI_ANOVA$ANOVA
PMI_ANOVA$`Mauchly's Test for Sphericity`
PMI_ANOVA$`Sphericity Corrections`
aovPMI <- aovEffectSize(ezObj = PMI_ANOVA, effectSize = "pes")
aovDispTable(aovPMI)

###### PLOTTING ######
## PMI ## 
# make control data-frame
control_PMI <- subset(PMIdata, PMIdata$DIAGNOSIS == 'HC')
control_PMI$TSTAT <- 0
control_PMI$PVALUE <- 1
control_PMI$DEFICIT <- 0
control_PMI$BL <- 0
# get deficit data for patients
plot_PMI <- merge(PMIdata, td_results, by = c('PPT','DOM','DIAGNOSIS'))
# include only relevant info
plot_PMI$PMI <- plot_PMI$PMI.x
plot_PMI <- plot_PMI[, c(1:10,13:14,22:24)]

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
  scale_shape_manual(values = c(1,18)) +
  facet_wrap(~DIAGNOSIS) +
  labs(x = 'Side', y = 'Lateral PMI (mm)') +
  theme_classic() + theme(legend.position = 'none', 
                          axis.title = element_text(size = 12),
                          axis.text = element_text(size = 10),
                          strip.text = element_text(size = 10) 
  ) -> pPMI
pPMI

## PLOT 4: average PMI across sides - combine with PLOT 2
PMIav_plot <- aggregate(PMI ~ PPT*DIAGNOSIS*AGE*ED*SITE, mean, data = plot_PMI)
PMIav_plot <- PMIav_plot[order(PMIav_plot$PPT),]
plot_PMI$DEFICITS <- as.numeric(plot_PMI$DEFICITS)
deficit <- aggregate(DEFICITS ~ PPT*DIAGNOSIS, max, data = plot_PMI)
deficit <- deficit[order(deficit$PPT),]
PMIav_plot$DEFICITS <- factor(deficit$DEFICITS)

ggplot(PMIav_plot, aes(x = DIAGNOSIS, y = PMI, colour = DIAGNOSIS, group = PPT, shape = DEFICITS)) + 
  geom_point(size = 2.5, position = position_dodge(.2)) +
  stat_summary(aes(y = PMI, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1, size = 4, group = 1) +
  scale_color_manual(values = c('grey45','grey45','grey45')) +
  scale_shape_manual(values = c(1,18)) +
  labs(x = 'Diagnosis', y = '') +
  theme_classic() + theme(legend.position = 'none', 
                          axis.title = element_text(size = 12),
                          axis.text = element_text(size = 10),
                          strip.text = element_text(size = 10) 
  ) -> avPMI
avPMI

PMIfig <- ggarrange(pPMI, avPMI,
                    ncol=2, nrow=1,
                    widths = c(1.5,1),
                    labels = c('a','b'),
                    hjust = -1)
PMIfig
ggsave('LATPMI-fig.png', plot = last_plot(), device = NULL, dpi = 300, 
       width = 8, height = 4, path = anaPath)


