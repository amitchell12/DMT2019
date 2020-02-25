library(readr)
library(ggplot2)
library(Rmisc)
library(gridExtra)
library(plyr)
library(tidyverse)
library(reshape2)
library(Hmisc)
library(ggpubr)

#set working directory to where data is
#on mac
#anaPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/analysis/lateral_reaching'
#dataPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/data'
#on pc
dataPath <- 'S:/groups/DMT/pilot-data'
anaPath <- 'S:/groups/DMT/pilot-data/analysis'
setwd(dataPath)

########### variable info ###########
visAngle <- function(size, distance){
  # this function calculates visual angle
  # size and distance must be in the same units
  Rad = 2*atan(size/(2*distance))
  Ang = Rad*(180/pi)
  return(Ang)
}
# screen information
x = 310
y = 175
pixels_perdeg = 43.76
deg_perpix = 1/pixels_perdeg
x_res = 1920
y_res = 1080
sd = 40

pixPer_mm_x = x_res/x;
pixPer_mm_y = y_res/y;
pixPer_mm = (pixPer_mm_x+pixPer_mm_y)/2;
mm_perPix = 1/pixPer_mm;

############ file info ############
#getting all datafiles and compiling (patient + control)
filenames <- dir(dataPath, recursive = TRUE, full.names = FALSE, pattern = '.csv')
# making results file
res <- read.csv(text = "subject_nr,targ_x,targ_y,land_x,land_y,target_onset,touch_offset,reach_duration,eye_move,void,task")

# name of csv file
# get key information from each csv file and then compile
### contine tomorrow/email Rob
for (file in filenames){
  if (isTRUE(substr(basename(file), 8, 8)=="8")){
    tmp <- read.csv(file)[, c(14:16,11:12,17,19,13,9,22)]
    tmp$task <- substr(basename(file), 25, 25)
    res <- rbind(res, tmp)
  }
}

# removing df (subject number = 300)
res <- res[res$subject_nr != 300, ]
res <- res[res$subject_nr != 500, ]
res <- res[res$land_x > -700, ]

# adding key details to data-frame
res$task <- factor(res$task) # task = factor
res$task = revalue(res$task, c('C'='closed', 'O'='open')) #renaming task
res$side <- factor(res$targ_x > 0, label = c('right')) #adding factor of side
res$ecc <- factor(cut(abs(res$targ_x), 3), labels = c('28', '33', '38')) #adding eccentricity
res$height <- factor(cut(res$targ_y, 3, labels = c('top', 'mid', 'bottom'))) #adding target height
#finally - target location
res$target <- paste0(res$ecc, res$height)
res$group <- factor(substr(res$subject_nr, 1, 1)) # adding group: 1 = HC, 2/4 = patient
#res$site <- factor(substr(res$subject_nr, 1, 1)) # 1.2 edinburgh, 3/4 norwich
# group '4' = patients from norwich, change to group = 2
#for (i in 1:length(res$group)){
#  if (isTRUE(res$group[i] == "8")){
#    res$group[i] = 2
#  }
#}
#for (i in 1:length(res$site)){
#  if (isTRUE(res$site[i] == "2")){
#    res$site[i] = 1
#  }
#}
res$group <- factor(res$group)
#res$site <- factor(res$site)

# eye move and void trials
# counting eye-move per participant
nEye_move <- aggregate(eye_move ~ subject_nr * group, sum, data = res)
totEye_move <- summarySE(nEye_move, measurevar = 'eye_move', groupvars = 'group')
nVoid <- aggregate(res$void, by=list(subject_nr = res$subject_nr), FUN=sum)

# calculating x and y error for each targ location
res$xerr_mm = (res$land_x - res$targ_x)*mm_perPix # in mm
res$yerr_mm = (res$land_y - res$targ_y)*mm_perPix
res$xerr_deg = visAngle(size= res$xerr_mm, distance= 400) # in deg
res$yerr_deg = visAngle(size= res$yerr_mm, distance= 400)

#absolute error in mm
res$AEmm = sqrt(res$xerr_mm^2 + res$yerr_mm^2)
res$AEdeg = sqrt(res$xerr_deg^2 + res$yerr_deg^2)

# removing and reorganising
res <- res[which(res$eye_move == 0 & res$void == 0), c(1,11:16,2:5,17:22,6:8)]

setwd(anaPath)
write.csv(res, "open-loop_data.csv", row.names = FALSE)

ggplot(res) + geom_point(aes(x = targ_x, y = targ_y), shape = 4, size = 3) +
  geom_point(aes(x = land_x, y = land_y, colour = ecc), shape = 1, size = 2) +
  facet_wrap(. ~subject_nr*task) -> allPP_plot

######### data aggregation + plotting ############
res_medians <- aggregate(AEdeg ~ ecc * side * task * subject_nr * group, median, data = res)
colnames(res_medians)[colnames(res_medians)=='AEdeg'] <- 'AEmed' #change name to be more logical
res_means <- aggregate(AEmed ~ task * side * subject_nr * group, mean, data = res_medians)
colnames(res_means)[colnames(res_means) == 'AEmed'] <- 'AEmean'
# save data
write.csv(res_medians, 'open-loop_medians.csv', row.names = FALSE)
# to calculate PMI need to cast by task....
PMIdata <- dcast(res_means, subject_nr+group+side ~ task) #different data-frame


# changing levels of PMI for plotting
res_means$side <- factor(res_means$side, levels = c('left', 'right'))
levels(res_means$side) <- c('Left', 'Right')
levels(res_means$group) <- c('Control', 'Patient') #changing group name from 1 = control, 2 = AD
levels(res_means$task) <- c('Peripheral', 'Free')
levels(res_means$site) <- c('UOE', 'UEA')
write.csv(res_means, 'lateral-reaching_means.csv', row.names = FALSE)

# mean plot 
ggplot(res_means, aes(x = task, y = AEmean)) +
  geom_point(shape = 1, size = 2) +
  geom_line(aes(group = subject_nr), size = 0.5, alpha = .5) + ylim(-.5,8) +
  labs(x = 'Side', y = 'Mean AE (deg)', element_text(size = 12)) +
  theme_bw() + theme(legend.position = 'none', text = element_text(size = 10),
                     strip.text.x = element_text(size = 8)) -> meansPlot

ggsave('allmeans_plot.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

######## not really needed rn
# summary
meanPMI_side <- summarySE(PMIdata, measurevar = 'PMI', groupvar = c('group', 'side'),
                       na.rm = TRUE)
meanPMI_all <- summarySE(PMIdata, measurevar = 'PMI', groupvar = c('group'),
                         na.rm = TRUE)
mean_periph_all <- summarySE(res_periph, measurevar = 'AEmean', groupvar = c('group'),
                             na.rm = TRUE)

# plot by eccentricity
# controls
meds_control <- res_medians[res_medians$group == 1 ,]
  
ggplot(meds_control, aes(x = ecc, y = AEmed, colour = side)) +
  geom_point(shape = 1, size = 2) +
  geom_line(aes(group = subject_nr), size = 0.5, alpha = .5) +
  facet_grid(cols = vars(side), rows = vars(task)) + ylim(-.5,10) +
  labs(title = 'Control', x = 'Eccentricity (deg)', 
       y = 'Mean AE (deg)', element_text(size = 12)) +
  scale_colour_manual(values = c('black', 'grey50')) +
  theme_bw() + theme(legend.position = 'none', text = element_text(size = 10),
                    strip.text.x = element_text(size = 10)) -> eccPlot


ggsave('control_ecc.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

# patients
meds_patient <- res_medians[res_medians$group == 2 ,]

ggplot(meds_patient, aes(x = ecc, y = AEmed, colour = side)) +
  geom_point(shape = 1, size = 2) +
  geom_line(aes(group = subject_nr), size = 0.5, alpha = .5) +
  facet_grid(cols = vars(side), rows = vars(task)) + ylim(-.5,12) +
  labs(title = 'Patient', x = 'Eccentricity (deg)', 
       y = 'Mean AE (deg)', element_text(size = 12)) +
  scale_colour_manual(values = c('black', 'grey50')) +
  theme_bw() + theme(legend.position = 'none', text = element_text(size = 10),
                     strip.text.x = element_text(size = 10)) -> eccPlot


ggsave('patient_ecc.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

###### directional error calc ######
dir_medians <- aggregate(xerr_deg ~ ecc * side * task * subject_nr * site * group, 
                         median, data = res)
colnames(dir_medians)[colnames(dir_medians)=='xerr_deg'] <- 'xerr_med' #change name to be more logical
dir_means <- aggregate(xerr_med ~ task * side * subject_nr * site * group, 
                       mean, data = dir_medians)
colnames(dir_means)[colnames(dir_means) == 'xerr_med'] <- 'xerr_mean'

# PMI for directional data (DMI - directional misreaching index)
dPMIdata <- dcast(dir_means, subject_nr+group+site+side ~ task) #different data-frame
dPMIdata$PMI <- dPMIdata$periph - dPMIdata$free
write.csv(dPMIdata, 'lateral-reaching_dPMI.csv', row.names = FALSE)

# plotting this
ggplot(dir_means, aes(x = task, y = xerr_mean, colour = site)) +
  geom_point(shape = 1, size = 2) +
  geom_line(aes(group = subject_nr), size = 0.5, alpha = .5) +
  facet_grid(cols = vars(side), rows = vars(group)) + ylim(-8,8) +
  labs(x = 'Side', y = 'Directional error (deg)', element_text(size = 12)) +
  scale_colour_manual(values = c('grey50', 'black')) +
  theme_bw() + theme(legend.position = 'none', text = element_text(size = 10),
                     strip.text.x = element_text(size = 8)) -> meansPlot

ggsave('directional_means_plot.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

# dPMI plot 
ggplot(dPMIdata, aes(x = side, y = PMI, colour = group), position = position_dodge(.2)) + 
  geom_point(shape = 1, size = 1.5, stroke = .8) +
  geom_line(aes(group = subject_nr), alpha = .5, size = .5) +
  scale_colour_manual(values = c('grey40', 'grey40')) +
  stat_summary(aes(y = PMI, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1, size = 3, group = 1) +
  ylim(-8,8) + labs(title = 'Lateral Reaching', x = 'Side', y = 'dPMI (deg)', 
                     element_text(size = 12)) +
  facet_wrap(~group) +
  theme_bw() + theme(legend.position = 'none', text = element_text(size = 10),
                     strip.text.x = element_text(size = 10)) -> dPMIplot

ggsave('dPMI_plot.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

## summary dPMI
meandPMI <- summarySE(dPMIdata, measurevar = 'PMI', groupvar = c('group', 'site', 'side'),
                          na.rm = TRUE)
meandPMI_all <- summarySE(dPMIdata, measurevar = 'PMI', groupvar = c('group', 'site', 'side'),
                         na.rm = TRUE)

##### outlier calculation, for controls only
controlData <- PMIdata[PMIdata$subject_nr < 200, ]

# median values for each side
tmp <- aggregate(controlData$PMI,  by=list(side = controlData$side), FUN=median)
names(tmp)[2] <- 'med'
controlData <- merge(tmp, controlData)

# calculating MAD for each (absolute value)
controlData$AD <- abs(controlData$PMI - controlData$med)
tmp <- aggregate(controlData$AD,  by=list(side = controlData$side), FUN=median)
names(tmp)[2] <- 'MAD'
controlData <- merge(controlData, tmp)

# adjusted z-score from these values
controlData$az <- (controlData$PMI - controlData$med)/(controlData$MAD * 1.4826)
controlData$z <- scale(controlData$PMI)

plot_name = 'adjustedZ.png'
AZplot <- ggplot(all_PMI, aes(x = side, y = az, colour = sub)) +
  geom_point(size = 3, position = position_dodge(.1)) + ylim(-3,6) + 
  stat_summary(aes(y = az, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', group = 1) + 
  labs(x = '', y = 'Adjusted z-score', element_text(size = 13)) +
  theme_bw() + theme(legend.position = "right", legend.title = element_blank())

ggsave(plot_name, plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, width = 5, height = 6.5, path = anaPath)


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

# plotting means
# offset
ggplot(res_offset_means, aes(x = side, y = time_touch_offset, colour = site)) +
  geom_point(shape = 1, size = 2) +
  geom_line(aes(group = subject_nr), size = 0.5, alpha = .5) +
  facet_grid(cols = vars(task), rows = vars(group)) + ylim(0, 1000) +
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
  facet_grid(cols = vars(task), rows = vars(group)) + ylim(0, 1000) +
  labs(title = 'Reach Duration', x = 'Side', 
       y = 'Reach duration (ms)', element_text(size = 12)) +
  scale_colour_manual(values = c('grey50', 'black')) +
  theme_bw() + theme(legend.position = 'none', text = element_text(size = 10),
                     strip.text.x = element_text(size = 10)) -> reachPlot

ggsave('reachDur_meansPlot.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, width = 5, height = 6.5, path = anaPath)

#means of both sides
res_reach_meansall <- aggregate(reach_duration~task * subject_nr * site * group, mean, 
                                data= res_reach_means) 
res_reach_meansall$task <- with(res_reach_meansall, factor(task, levels = rev(levels(task))))

ggplot(res_reach_meansall, aes(x = task, y = reach_duration, colour = site)) +
  geom_point(shape = 1, size = 2) +
  geom_line(aes(group = subject_nr), size = 0.5, alpha = .5) +
  facet_wrap(~group) + ylim(0, 1000) +
  labs(title = 'Lateral reaching', x = 'Task', 
       y = 'Reach duration (ms)', element_text(size = 12)) +
  scale_colour_manual(values = c('grey50', 'black')) +
  theme_bw() + theme(legend.position = 'none', text = element_text(size = 10),
                     strip.text.x = element_text(size = 10)) -> reachPlot

ggsave('reachDur_meansPlot.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)


# correlating peripheral reach duration with PMI
# cast 
rt_offset <- dcast(res_offset_means, subject_nr+group+site+side ~ task)
rt_reach <- dcast(res_reach_means, subject_nr+group+site+side ~ task) #different data-frame
#correlations
corrData <- data.frame(PMIdata$PMI)
colnames(corrData)[colnames(corrData) == 'PMIdata.PMI'] <- 'PMI' #renaming
corrData$pAE <- PMIdata$periph
corrData$reachRT <- rt_reach$periph
corrData$offsetRT <- rt_offset$periph
corrData$group <- PMIdata$group
corrData$side <- PMIdata$side

# reach dur plot
ggscatter(corrData, x = "pAE", y = "reachRT", add = 'reg.line', conf.int = TRUE,
          cor.coef = TRUE, cor.method = 'pearson') + 
  facet_grid(cols = vars(side), rows = vars(group))
# touch offset plot
ggscatter(corrData, x = "pAE", y = "offsetRT", add = 'reg.line', conf.int = TRUE,
          cor.coef = TRUE, cor.method = 'pearson') + 
  facet_grid(cols = vars(side), rows = vars(group))


########### next steps: comparing patients to controls




### isolating P12
for (i in 1:length(res_means$task)){ 
  if (isTRUE(res_means$subject_nr[i] == '212')){ 
    res_means$group <- as.numeric(res_means$group)
    res_means$group[i] = 3}
  }
res_means$group <- factor(res_means$group)

## plotting
ggplot(res_means, aes(x = side, y = AEmean, colour = group)) +
  geom_point(shape = 1, size = 2) +
  geom_line(aes(group = subject_nr), size = 0.5, alpha = .5) +
  facet_grid(cols = vars(task), rows = vars(group)) + ylim(-.5,8) +
  labs(x = 'Side', y = 'Mean AE (deg)', element_text(size = 12)) +
  scale_colour_manual(values = c('black', 'grey50', 'red')) +
  theme_bw() + theme(legend.position = 'none', text = element_text(size = 10),
                     strip.text.x = element_text(size = 8)) -> meansPlot


## PMI
