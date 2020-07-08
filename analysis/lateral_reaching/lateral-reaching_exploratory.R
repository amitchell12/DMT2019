library(readr)
library(ggplot2)
library(Rmisc)
library(gridExtra)
library(plyr)
library(tidyverse)
library(reshape2)
library(Hmisc)
library(ggpubr)

###### GETTING DATA #######
#on mac
anaPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/analysis/lateral_reaching'
# on desktop mac
#anaPath <- '/Users/Alex/Documents/DMT/analysis/lateral_reaching'
#dataPath <- '/Users/Alex/Documents/DMT/data'
#on pc
#dataPath <- 'S:/groups/DMT/data'
#anaPath <- 'S:/groups/DMT/analysis/lateral_reaching'
setwd(anaPath)

res <- read.csv('lateral-reaching_compiled.csv')
## getting dominant + non-dominant sides, for analysis
names(res)[4] <- 'HAND'
res$HAND <- factor(res$HAND, labels = c('left','right'))
# if hand = side, dominant; else non dominant
res$DOM <- as.numeric(res$SIDE == res$HAND) #1 = dominant, 0 = non-dominant
res$DOM <- factor(res$DOM, labels= c('ND','D'))
# change order so dominance up-front
res <- res[, c(1:8,27,9:26)]

# changing levels to be more informative
res$GRP <- factor(res$GRP, labels = c('Control','Patient'))
levels(res$VIEW) <- c('Free', 'Peripheral')
res$SITE <- factor(res$SITE, labels = c('UOE','UEA'))
res$DIAGNOSIS <- factor(res$DIAGNOSIS)

## find outliers and remove ##
xclude <- read.csv('lateraloutliers.csv')
res <- res[!(res$PPT %in% xclude$PPT), ]

###### DIRECTIONAL ERROR: PMI ######
dir_medians <- aggregate(xerr_mm ~ PPT * VIEW * SIDE * POSITION * SITE * GRP * DIAGNOSIS, 
                         median, data = res)
colnames(dir_medians)[colnames(dir_medians)=='xerr_mm'] <- 'xerr_med' #change name to be more logical
dir_means <- aggregate(xerr_med ~ PPT* VIEW * SIDE * SITE * GRP * DIAGNOSIS, 
                       mean, data = dir_medians)
colnames(dir_means)[colnames(dir_means) == 'xerr_med'] <- 'xerr_mean'

# PMI for directional data (DMI - directional misreaching index)
dPMIdata <- dcast(dir_means, PPT+SIDE+DIAGNOSIS+SITE ~ VIEW) #different data-frame
dPMIdata$PMI <- dPMIdata$Peripheral - dPMIdata$Free
write.csv(dPMIdata, 'lateral-reaching_dirPMI.csv', row.names = FALSE)

# plotting this
ggplot(dir_means, aes(x = VIEW, y = xerr_mean, colour = DIAGNOSIS)) +
  geom_point(shape = 1, size = 2) +
  geom_line(aes(group = PPT), size = 0.5, alpha = .5) +
  facet_grid(cols = vars(SIDE), rows = vars(DIAGNOSIS)) + 
  labs(x = '', y = 'Directional error (mm)', element_text(size = 12)) +
  theme_bw() + theme(legend.position = 'none', text = element_text(size = 10),
                     strip.text.x = element_text(size = 8)) 

# dPMI plot 
ggplot(dPMIdata, aes(x = SIDE, y = PMI, colour = DIAGNOSIS)) + 
  geom_point(shape = 1, size = 1.5, stroke = .8) +
  geom_line(aes(group = PPT), alpha = .5, size = .5) +
  stat_summary(aes(y = PMI, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1, size = 3, group = 1) +
  labs(title = 'Lateral Reaching', x = 'Side', y = 'dPMI (mm)', 
                    element_text(size = 12)) +
  facet_wrap(~DIAGNOSIS) +
  theme_bw() + theme(legend.position = 'none', text = element_text(size = 10),
                     strip.text.x = element_text(size = 10))

ggsave('dPMI_plot.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

## summary dPMI
mean_dPMI <- summarySE(dPMIdata, measurevar = 'PMI', groupvar = c('DIAGNOSIS', 'SIDE'),
                       na.rm = TRUE)
mean_dPMI_all <- summarySE(dPMIdata, measurevar = 'PMI', groupvar = c('DIAGNOSIS', 'SIDE'),
                           na.rm = TRUE)

##### DIRECTIONAL ERROR: ECCENTRICITY #####


### analysis of response times
# do in the same manner as with absolute error
# medians
res_offset_medians <- aggregate(
  time_touch_offset ~ ecc * side * task * subject_nr * site * group, median, data = res)
res_reach_medians <- aggregate(
  reach_duration ~ ecc * side * task * subject_nr * site * group, median, data = res)
# means of medians
res_offset_means <- aggregate(
  time_touch_offset ~ task * side * subject_nr * site * group, mean, data = res_offset_medians)
res_reach_means <- aggregate(
  reach_duration ~ task * side * subject_nr * site * group, mean, data = res_reach_medians)
levels(res_offset_means$group) <- c('Control', 'Patient')
levels(res_reach_means$group) <- c('Control', 'Patient')
levels(res_offset_means$site) <- c('UOE', 'UEA')
levels(res_reach_means$site) <- c('UOE', 'UEA')
#adding diagnosis to this info
res_offset_means <- merge(demo, res_offset_means)
res_reach_means <- merge(demo, res_reach_means)


# plotting means
# offset
ggplot(res_offset_means, aes(x = side, y = time_touch_offset, colour = site)) +
  geom_point(shape = 1, size = 2) +
  geom_line(aes(group = subject_nr), size = 0.5, alpha = .5) +
  facet_grid(cols = vars(task), rows = vars(diagnosis)) + ylim(0, 1000) +
  labs(title = 'Touch Offset', x = 'Side', 
       y = 'Touch offset RT (ms)', element_text(size = 12)) +
  scale_colour_manual(values = c('grey50', 'black')) +
  theme_bw() + theme(legend.position = 'none', text = element_text(size = 10),
                     strip.text.x = element_text(size = 10)) -> touchoffsetPlot

ggsave('touchoffset_meansPlot.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, width = 5, height = 6.5, path = anaPath)
  
#reach
ggplot(res_reach_means, aes(x = side, y = reach_duration, colour = site)) +
  geom_point(shape = 1, size = 2) +
  geom_line(aes(group = subject_nr), size = 0.5, alpha = .5) +
  facet_grid(cols = vars(task), rows = vars(diagnosis)) + ylim(0, 1000) +
  labs(title = 'Reach Duration', x = 'Side', 
       y = 'Reach duration (ms)', element_text(size = 12)) +
  scale_colour_manual(values = c('grey50', 'black')) +
  theme_bw() + theme(legend.position = 'none', text = element_text(size = 10),
                     strip.text.x = element_text(size = 10)) -> reachPlot

ggsave('reachDur.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, width = 5, height = 6.5, path = anaPath)

#means of both sides
res_reach_meansall <- aggregate(reach_duration~task * subject_nr * site * diagnosis, mean, 
                                data= res_reach_means) 
res_reach_meansall$task <- with(res_reach_meansall, factor(task, levels = rev(levels(task))))

ggplot(res_reach_meansall, aes(x = task, y = reach_duration, colour = site)) +
  geom_point(shape = 1, size = 2) +
  geom_line(aes(group = subject_nr), size = 0.5, alpha = .5) +
  facet_wrap(~diagnosis) + ylim(0, 1000) +
  labs(title = 'Lateral reaching', x = 'Task', 
       y = 'Reach duration (ms)', element_text(size = 12)) +
  scale_colour_manual(values = c('grey50', 'black')) +
  theme_bw() + theme(legend.position = 'none', text = element_text(size = 10),
                     strip.text.x = element_text(size = 10)) -> reachPlot

ggsave('reachDur_means.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)


# correlating peripheral reach duration with PMI
# cast 
rt_offset <- dcast(res_offset_means, subject_nr+diagnosis+site+side ~ task)
rt_reach <- dcast(res_reach_means, subject_nr+diagnosis+site+side ~ task) #different data-frame
#correlations
corrData <- data.frame(PMIdata$PMI)
colnames(corrData)[colnames(corrData) == 'PMIdata.PMI'] <- 'PMI' #renaming
corrData$pAE <- PMIdata$periph
corrData$reachRT <- rt_reach$periph
corrData$offsetRT <- rt_offset$periph
corrData$diagnosis <- PMIdata$diagnosis
corrData$side <- PMIdata$side

# reach dur plot
ggscatter(corrData, x = "pAE", y = "reachRT", add = 'reg.line', conf.int = TRUE,
          cor.coef = TRUE, cor.method = 'pearson') + 
  facet_grid(cols = vars(side), rows = vars(diagnosis))
# touch offset plot
ggscatter(corrData, x = "pAE", y = "offsetRT", add = 'reg.line', conf.int = TRUE,
          cor.coef = TRUE, cor.method = 'pearson') + 
  facet_grid(cols = vars(side), rows = vars(diagnosis))


########### next steps: comparing patients to controls
