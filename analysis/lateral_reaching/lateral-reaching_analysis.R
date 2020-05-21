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

#set working directory to where data is
#on mac
#dataPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/data'
#anaPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/analysis/lateral_reaching'
# on desktop mac
anaPath <- '/Users/Alex/Documents/DMT/analysis/lateral_reaching'
dataPath <- '/Users/Alex/Documents/DMT/data'
#on pc
#dataPath <- 'S:/groups/DMT/data'
#anaPath <- 'S:/groups/DMT/analysis/lateral_reaching'
setwd(anaPath)

# load data file
res <- read.csv('lateral-reaching_compiled.csv')

######### step 1, calculating PMI, all data ############
res_medians <- aggregate(
  AEdeg ~ PPT * SIDE * VIEW * POSITION * SITE * GRP * DIAGNOSIS * AGE * ED, 
  median, data = res)
colnames(res_medians)[colnames(res_medians)=='AEdeg'] <- 'AEmed' #change name to be more logical

# changing levels to be more informative
res_medians$SIDE <- factor(res_medians$SIDE, levels = c('left', 'right'))
res_medians$GRP <- factor(res_medians$GRP)
levels(res_medians$GRP) <- c('Control', 'Patient') #changing group name from 1 = control, 2 = AD
levels(res_medians$VIEW) <- c('Free', 'Peripheral')
res_medians$SITE <- factor(res_medians$SITE)
levels(res_medians$SITE) <- c('UOE', 'UEA')
res_medians$DIAGNOSIS <- factor(res_medians$DIAGNOSIS)
res_medians <- res_medians[order(res_medians$PPT), ] 

res_means <- aggregate(AEmed ~ PPT * SIDE * VIEW * SITE * GRP * DIAGNOSIS * AGE * ED, 
                       mean, data = res_medians)
res_means <- res_means[order(res_means$PPT), ] 
colnames(res_means)[colnames(res_means) == 'AEmed'] <- 'AEmean'
# save data
write.csv(res_medians, 'lateral-reaching_medians.csv', row.names = FALSE)
write.csv(res_means, 'lateral-reaching_means.csv', row.names = FALSE)

# to calculate PMI need to cast by task....
PMIdata <- dcast(res_means, PPT+GRP+SITE+SIDE+DIAGNOSIS+AGE+ED ~ VIEW) #different data-frame
PMIdata$PMI <- PMIdata$Peripheral - PMIdata$Free
write.csv(PMIdata, 'lateral-reaching_PMI.csv', row.names = FALSE)


# mean plot 
ggplot(res_means, aes(x = SIDE, y = AEmean, colour = SITE)) +
  geom_point(shape = 1, size = 2) +
  geom_line(aes(group = PPT), size = 0.5, alpha = .5) +
  facet_grid(cols = vars(VIEW), rows = vars(DIAGNOSIS)) + ylim(-.5,8) +
  labs(x = 'Side', y = 'Mean AE (deg)', element_text(size = 12)) +
  scale_colour_manual(values = c('grey50', 'black')) +
  theme_bw() + theme(legend.position = 'none', text = element_text(size = 10),
                     strip.text.x = element_text(size = 8)) -> meansPlot

ggsave('allmeans_plot.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

# PMI plot 
ggplot(PMIdata, aes(x = SIDE, y = PMI, colour = SITE), position = position_dodge(.2)) + 
  geom_point(shape = 1, size = 4) +
  geom_line(aes(group = PPT), alpha = .5, size = .8) +
  scale_colour_manual(values = c('grey50', 'black')) +
  stat_summary(aes(y = PMI, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1, size = 5, group = 1) +
  ylim(-.5,10) + labs(title = '', x = 'Side', y = 'Reaching error (deg)', 
                     element_text(size = 14)) +
  facet_wrap(~DIAGNOSIS) +
  theme_bw() + theme(legend.position = 'none', text = element_text(size = 14),
                     strip.text.x = element_text(size = 12)) -> PMIplot

ggsave('lateralPMI.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, width = 7, height = 4, path = anaPath)

# summary
meanPMI_side <- summarySE(PMIdata, measurevar = 'PMI', groupvar = c('DIAGNOSIS', 'SIDE'),
                       na.rm = TRUE)
meanPMI_all <- summarySE(PMIdata, measurevar = 'PMI', groupvar = c('DIAGNOSIS'),
                         na.rm = TRUE)

# plot by eccentricity
# controls
meds_control <- res_medians[res_medians$GRP == 'Control' ,]
  
ggplot(meds_control, aes(x = POSITION, y = AEmed, colour = SIDE)) +
  geom_point(shape = 1, size = 2) +
  geom_line(aes(group = PPT), size = 0.5, alpha = .5) +
  facet_grid(cols = vars(SIDE), rows = vars(VIEW)) + ylim(-.5,10) +
  labs(title = 'Control', x = 'Eccentricity (deg)', 
       y = 'Mean AE (deg)', element_text(size = 12)) +
  scale_colour_manual(values = c('black', 'grey50')) +
  theme_bw() + theme(legend.position = 'none', text = element_text(size = 10),
                    strip.text.x = element_text(size = 10)) 


ggsave('control_ecc.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

# patients - MCI
meds_MCI <- res_medians[res_medians$DIAGNOSIS == 'MCI' ,]

ggplot(meds_MCI, aes(x = POSITION, y = AEmed, colour = SIDE)) +
  geom_point(shape = 1, size = 2) +
  geom_line(aes(group = PPT), size = 0.5, alpha = .5) +
  facet_grid(cols = vars(SIDE), rows = vars(VIEW)) + ylim(-.5,12) +
  labs(title = 'MCI', x = 'Eccentricity (deg)', 
       y = 'Mean AE (deg)', element_text(size = 12)) +
  scale_colour_manual(values = c('black', 'grey50')) +
  theme_bw() + theme(legend.position = 'none', text = element_text(size = 10),
                     strip.text.x = element_text(size = 10)) -> MCIecc

ggsave('MCI_ecc.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

# patients - AD
meds_AD <- res_medians[res_medians$DIAGNOSIS == 'AD' ,]

ggplot(meds_AD, aes(x = POSITION, y = AEmed, colour = SIDE)) +
  geom_point(shape = 1, size = 2) +
  geom_line(aes(group = PPT), size = 0.5, alpha = .5) +
  facet_grid(cols = vars(SIDE), rows = vars(VIEW)) + ylim(-.5,12) +
  labs(title = 'Alzheimers', x = 'Eccentricity (deg)', 
       y = 'Mean AE (deg)', element_text(size = 12)) +
  scale_colour_manual(values = c('black', 'grey50')) +
  theme_bw() + theme(legend.position = 'none', text = element_text(size = 10),
                     strip.text.x = element_text(size = 10)) 

ggsave('AD_ecc.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

# PMI collapsed across side
PMIav <- aggregate(PMI ~ PPT * SITE * GRP * DIAGNOSIS * AGE * ED, mean, data = PMIdata)
PMIav <- PMIav[order(PMIav$PPT), ]

######## step 2: single case stats, all controls #########
# correlate PMI with possible co-variates, age + education
agecov <- cor.test(PMIdata$AGE, PMIdata$PMI, method = 'pearson')
ggscatter(PMIdata, x = "AGE", y = "PMI", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Age", ylab = "PMI (deg")

# correlate PMI with years of education
PMIdata$ED <- as.numeric(PMIdata$ED)
edcov <- cor.test(PMIdata$ED, PMIdata$PMI, method = 'pearson')
ggscatter(PMIdata, x = "ED", y = "PMI", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Years of Education", ylab = "PMI (deg")

# create data-frames (controls and patients) with key information
td_summary <- summarySE(data = PMIdata, measurevar = 'PMI', 
                        groupvars = c('DIAGNOSIS','SIDE'), na.rm = TRUE)
# isolating control data for analysis
td_control <- td_summary[td_summary$DIAGNOSIS == 'HC', ]

#patient data for test of deficit
td_patient <- PMIdata[, c(1,4,5,10)]
td_patient <- td_patient[td_patient$DIAGNOSIS != 'HC' ,] #removing controls
td_patient <- dcast(td_patient, PPT+DIAGNOSIS ~ SIDE)
# NA values  = -1, so t-value is negative and can remove later
td_patient[is.na(td_patient$right), "right"] <- -1 

# time for test of deficit! Calling on 'singcar' package developed by Jonathan Rittmo
# using Crawford's (1998) test of deficit
td_right <- read.csv(text = 'PMI,TSTAT,PVALUE,SIDE,PPT,DIAGNOSIS')
td_left <- read.csv(text = 'PMI,TSTAT,PVALUE,SIDE,PPT,DIAGNOSIS')
for (l in 1:length(td_patient$PPT)){
  #left data first
  leftres <- TD(td_patient$left[l], td_control$PMI[1], td_control$sd[1], 24, 
            alternative = 'greater', na.rm = FALSE)
  ltmp <- data.frame(td_patient$left[l], leftres$statistic, leftres$p.value, 'left') 
  ltmp$PPT <- td_patient$PPT[l]
  ltmp$DIAGNOSIS <- td_patient$DIAGNOSIS[l]
  td_left <- rbind(td_left, ltmp)
  #then right data
  rightres <- TD(td_patient$right[l], td_control$PMI[2], td_control$sd[2], 24, 
                alternative = 'greater', na.rm = FALSE)
  rtmp <- data.frame(td_patient$right[l], rightres$statistic, rightres$p.value, 'right') 
  rtmp$PPT <- td_patient$PPT[l]
  rtmp$DIAGNOSIS <- td_patient$DIAGNOSIS[l]
  td_right <- rbind(td_right, rtmp)
}

# merging and renaming data-frames
td_results <- read.csv(text = 'PMI,TSTAT,PVALUE,SIDE,PPT,DIAGNOSIS')
# changing names of td data-frames to match td-res
names(td_left) <- names(td_results)
names(td_right) <- names(td_results)
td_results <- rbind(td_results, td_left, td_right)
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

# seperating data-frame into different groups
MCI <- td_results[td_results$DIAGNOSIS == 'MCI' ,]
AD <- td_results[td_results$DIAGNOSIS == 'AD' ,]
left <- td_results[td_results$SIDE == 'left' ,]
right <- td_results[td_results$SIDE == 'right' ,]

## refine plot later - need to read Crawford reporting
## plotting p-values
ggplot(td_results, aes(x = SIDE, y = PVALUE, group = PPT, colour = DIAGNOSIS)) +
  geom_point(size = 2, position = position_dodge(.2)) +
  geom_line(aes(group = PPT), size = 0.5, alpha = .5, position = position_dodge(.2)) + 
  scale_color_manual(values = c('black', 'grey50')) +
  geom_hline(yintercept = 0.025, color = 'red') +
  geom_hline(yintercept = 0.05, color = 'blue') +
  labs(title = 'Single case stats', x = 'Side', y = 'Probability of deficit', 
       element_text(size = 10)) +
  theme_bw()

ggsave('SCS_probability-scatter.png', plot = last_plot(), device = NULL, dpi = 300, 
       width = 6, height = 7, path = anaPath)

### binomial test ###
# calculating the likelihood that inflated PMI occurs at above chance in each patient group
# seperating out MCI and AD data frames into left and right
MCIleft <- MCI[MCI$SIDE == 'left' ,]
MCIright <- MCI[MCI$SIDE == 'right' ,]
ADleft <- AD[AD$SIDE == 'left' ,]
ADright <- AD[AD$SIDE == 'right' ,]

pval <- .05

## first do test of deficit
# binomial MCI
binMCI_left <- binom.test(sum(MCIleft$DEFICIT), length(MCIleft$PPT), pval, alternative = 'greater')
binMCI_right <- binom.test(sum(MCIright$DEFICIT), length(MCIright$PPT), pval, alternative = 'greater')
# binomial AD
binAD_left <- binom.test(sum(ADleft$DEFICIT), length(ADleft$PPT), pval, alternative = 'greater')
binAD_right <- binom.test(sum(ADright$DEFICIT), length(ADright$PPT), pval, alternative = 'greater')
# binomial all
bin_left <- binom.test(sum(left$DEFICIT), length(left$PPT), pval, alternative = 'greater')
bin_right <- binom.test(sum(right$DEFICIT), length(right$PPT), pval, alternative = 'greater')

## no cases that are borderline and not full deficit ( > 0.025, yet < 0.05), 
# so no need to run borderline analysis here

######## step 3: ANOVA, all controls #########

######## step 4: outlier removal, filtered PMI #########
controlData <- PMIdata[PMIdata$PPT < 200, ]

# median values for each side
tmp <- aggregate(PMI ~ GRP, median, data = controlData)
names(tmp)[2] <- 'med'
controlData <- merge(tmp, controlData)

# calculating MAD for each (absolute value)
controlData$AD <- abs(controlData$PMI - controlData$med)
tmp <- aggregate(AD ~ GRP, median, data=controlData)
names(tmp)[2] <- 'MAD'
controlData <- merge(controlData, tmp)

# adjusted z-score from these values
controlData$az <- (controlData$PMI - controlData$med)/(controlData$MAD * 1.4826)
controlData$z <- scale(controlData$PMI)
controlData$PPT <- factor(controlData$PPT)

plot_name = 'adjustedZ.png'
ggplot(controlData, aes(x = GRP, y = az, colour = PPT)) +
  geom_point(size = 3, position = position_dodge(.1)) +  
  stat_summary(aes(y = az, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', group = 1) + 
  labs(x = '', y = 'Adjusted z-score', element_text(size = 13)) +
  theme_bw() + theme(legend.position = "right", legend.title = element_blank())

ggsave(plot_name, plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, width = 5, height = 6.5, path = anaPath)

### outlier removal ###
# find controls with az > 2.5 - need to remove entire control, not just side
XCLUDE <- controlData[controlData$az > 2.5, ]
write.csv(XCLUDE, 'lateraloutliers.csv', row.names = FALSE)
# creating data-frame with control data removed
PMIfilt <- PMIdata[!(PMIdata$PPT %in% XCLUDE$PPT), ]
write.csv(PMIfilt, 'lateralPMI-filtered.csv', row.names = FALSE)

##### need to change the ordering of plot data-frames, worry about this for publication
### PLOT FILTERED PMI DATA
ggplot(PMIfilt, aes(x = SIDE, y = PMI, colour = SITE), position = position_dodge(.2)) + 
  geom_point(shape = 1, size = 4) +
  geom_line(aes(group = PPT), alpha = .5, size = .8) +
  scale_colour_manual(values = c('grey50', 'black')) +
  stat_summary(aes(y = PMI, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1, size = 5, group = 1) +
  ylim(-.5,10) + labs(title = 'Lateral Reaching', x = 'Side', y = 'Reaching error (deg)', 
                      element_text(size = 14)) +
  facet_wrap(~DIAGNOSIS) +
  theme_bw() + theme(legend.position = 'none', text = element_text(size = 14),
                     strip.text.x = element_text(size = 12)) 

ggsave('lateralPMI-filtered.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, width = 7, height = 4, path = anaPath)

## averaging PMI data across sides
meanFPMI <- summarySE(PMIfilter, measurevar = 'PMI', groupvar = c('DIAGNOSIS', 'SIDE'),
                      na.rm = TRUE)
meanFPMI_all <- summarySE(PMIfilter, measurevar = 'PMI', groupvar = c('DIAGNOSIS'),
                          na.rm = TRUE)

#average across side
PMIfilt_av <- aggregate(PMI ~ PPT * SITE * DIAGNOSIS, mean, data = PMIfilt)
jitter <- position_jitter(width = 0.1, height = 0.1)

ggplot(PMIfilt_av, aes(x = DIAGNOSIS, y = PMI, colour = SITE)) + 
  geom_point(position = jitter, shape = 21, size = 3) +
  scale_colour_manual(values = c('grey40', 'black')) +
  stat_summary(aes(y = PMI, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1, size = 4, group = 1) +
  ylim(-.5,7) + labs(title = '', x = '', y = 'Reaching error (deg)', 
                     element_text(size = 8)) +
  theme_bw() + theme(legend.position = 'none', text = element_text(size = 10),
                     strip.text.x = element_text(size = 8))

ggsave('lateralPMI-filtered-av.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, width = 3, height = 3, path = anaPath)

## ANOVA ## 


###### step xx DIRECTIONAL ERROR ######
dir_medians <- aggregate(xerr_deg ~ POSITION * SIDE * VIEW * PPT * SITE * GRP * DIAGNOSIS, 
                         median, data = res)
colnames(dir_medians)[colnames(dir_medians)=='xerr_deg'] <- 'xerr_med' #change name to be more logical
dir_means <- aggregate(xerr_med ~ VIEW * SIDE * PPT * SITE * GRP * DIAGNOSIS, 
                       mean, data = dir_medians)
colnames(dir_means)[colnames(dir_means) == 'xerr_med'] <- 'xerr_mean'

# PMI for directional data (DMI - directional misreaching index)
dPMIdata <- dcast(dir_means, PPT+GRP+DIAGNOSIS+SITE+SIDE ~ VIEW) #different data-frame
dPMIdata$PMI <- dPMIdata$periph - dPMIdata$free
write.csv(dPMIdata, 'lateral-reaching_dPMI.csv', row.names = FALSE)

# plotting this
ggplot(dir_means, aes(x = VIEW, y = xerr_mean, colour = SITE)) +
  geom_point(shape = 1, size = 2) +
  geom_line(aes(group = PPT), size = 0.5, alpha = .5) +
  facet_grid(cols = vars(SIDE), rows = vars(DIAGNOSIS)) + ylim(-8,8) +
  labs(x = 'Side', y = 'Directional error (deg)', element_text(size = 12)) +
  scale_colour_manual(values = c('grey50', 'black')) +
  theme_bw() + theme(legend.position = 'none', text = element_text(size = 10),
                     strip.text.x = element_text(size = 8)) -> meansPlot

ggsave('directional_means_plot.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

# dPMI plot 
ggplot(dPMIdata, aes(x = SIDE, y = PMI, colour = SITE), position = position_dodge(.2)) + 
  geom_point(shape = 1, size = 1.5, stroke = .8) +
  geom_line(aes(group = PPT), alpha = .5, size = .5) +
  scale_colour_manual(values = c('grey40', 'black')) +
  stat_summary(aes(y = PMI, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1, size = 3, group = 1) +
  ylim(-8,8) + labs(title = 'Lateral Reaching', x = 'Side', y = 'dPMI (deg)', 
                     element_text(size = 12)) +
  facet_wrap(~DIAGNOSIS) +
  theme_bw() + theme(legend.position = 'none', text = element_text(size = 10),
                     strip.text.x = element_text(size = 10)) -> dPMIplot

ggsave('dPMI_plot.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

## summary dPMI
mean_dPMI <- summarySE(dPMIdata, measurevar = 'PMI', groupvar = c('DIAGNOSIS', 'SIDE'),
                          na.rm = TRUE)
mean_dPMI_all <- summarySE(dPMIdata, measurevar = 'PMI', groupvar = c('DIAGNOSIS', 'SIDE'),
                         na.rm = TRUE)
