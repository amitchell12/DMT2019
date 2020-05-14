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
#anaPath <- '/Users/Alex/Documents/DMT/analysis/lateral_reaching'
#dataPath <- '/Users/Alex/Documents/DMT/data'
#on pc
dataPath <- 'S:/groups/DMT/data'
anaPath <- 'S:/groups/DMT/analysis/lateral_reaching'
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
                     strip.text.x = element_text(size = 10))

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

######## step 3: ANOVA, all controls #########

######## step 4: outlier removal, filtered PMI #########
controlData <- PMIav[PMIav$PPT < 200, ]

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
AZplot <- ggplot(controlData, aes(x = GRP, y = az, colour = PPT)) +
  geom_point(size = 3, position = position_dodge(.1)) + ylim(-2,4) + 
  stat_summary(aes(y = az, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', group = 1) + 
  labs(x = '', y = 'Adjusted z-score', element_text(size = 13)) +
  theme_bw() + theme(legend.position = "right", legend.title = element_blank())

ggsave(plot_name, plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, width = 5, height = 6.5, path = anaPath)

####### outlier removal #######
# find controls with az > 2.5 and flag-up in PMI data-set - need to remove entire control, not just side
for (l in 1:length(controlData$PPT)){
  if (isTRUE(controlData$az[l] > '2.5')){
    controlData$outlier[l] = 1
  }
  else 
    controlData$outlier[l] = 0
}
#reorder control data
controlData <- controlData[order(controlData$PPT), ]

# adding to PMI data set, then removing
##### reached here- try and merge outlier to PMI data, need a fresher brain for this... 
PMIfilter <- merge(controlData, PMIdata, by = 'PPT')
PMIfilter <- PMIfilter[, c(1:8,14)]
# save with outlier information
write.csv(PMIfilter, 'lateral-outliers.csv', row.names = FALSE)

## removing outliers
for (l in 49:length(PMIfilter$outlier)){
  PMIfilter$outlier[l] = 0
}
## this removes the data points, but not the controls
PMIfilter <- PMIfilter[PMIfilter$outlier < 1, ]

### PLOT FILTERED PMI DATA
ggplot(PMIfilter, aes(x = SIDE, y = PMI, colour = SITE), position = position_dodge(.2)) + 
  geom_point(shape = 1, size = 4) +
  geom_line(aes(group = PPT), alpha = .5, size = .8) +
  scale_colour_manual(values = c('grey50', 'black')) +
  stat_summary(aes(y = PMI, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1, size = 5, group = 1) +
  ylim(-.5,10) + labs(title = 'Lateral Reaching', x = 'Side', y = 'Reaching error (deg)', 
                      element_text(size = 14)) +
  facet_wrap(~DIAGNOSIS) +
  theme_bw() + theme(legend.position = 'none', text = element_text(size = 14),
                     strip.text.x = element_text(size = 12)) -> PMIf_plot

ggsave('lateralPMI-filtered.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, width = 7, height = 4, path = anaPath)

meanFPMI <- summarySE(PMIfilter, measurevar = 'PMI', groupvar = c('DIAGNOSIS', 'SIDE'),
                      na.rm = TRUE)
meanFPMI_all <- summarySE(PMIfilter, measurevar = 'PMI', groupvar = c('DIAGNOSIS'),
                          na.rm = TRUE)

#average across side
PMIfilter_av <- aggregate(PMI ~ PPT * SITE * DIAGNOSIS, mean, data = PMIfilter)
jitter <- position_jitter(width = 0.1, height = 0.1)

ggplot(PMIfilter_av, aes(x = DIAGNOSIS, y = PMI, colour = SITE)) + 
  geom_point(position = jitter, shape = 21, size = 3) +
  scale_colour_manual(values = c('grey40', 'black')) +
  stat_summary(aes(y = PMI, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1, size = 4, group = 1) +
  ylim(-.5,7) + labs(title = '', x = '', y = 'Reaching error (deg)', 
                     element_text(size = 8)) +
  theme_bw() + theme(legend.position = 'none', text = element_text(size = 10),
                     strip.text.x = element_text(size = 8)) -> PMIf_plot

ggsave('lateralPMI-filtered-av.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, width = 3, height = 3, path = anaPath)

## ANOVA ## 

###### step 3 SINGLE CASE STATS ######
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
