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
dataPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/data'
anaPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/analysis/lateral_reaching'
# on desktop mac
#anaPath <- '/Users/Alex/Documents/DMT/analysis/lateral_reaching'
#dataPath <- '/Users/Alex/Documents/DMT/data'
#on pc
#dataPath <- 'S:/groups/DMT/data'
#anaPath <- 'S:/groups/DMT/analysis/lateral_reaching'
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
  if (isTRUE(substr(basename(file), 12, 12)=="p")){
    tmp <- read.csv(file)[, c(14:16,11:12,17,19,13,9,22)]
    tmp$task <- substr(basename(file), 28, 28)
    res <- rbind(res, tmp)
  }
}

# removing df (subject number = 300, DF and ==500 PCA, Alicia)
res <- res[res$subject_nr != 300, ]
res <- res[res$subject_nr != 500, ]

# adding key details to data-frame
res$task <- factor(res$task) # task = factor
res$task = revalue(res$task, c('i'='periph', 'r'='free')) #renaming task
res$side <- factor(res$targ_x < 0, label = c('right', 'left')) #adding factor of side
res$ecc <- factor(cut(abs(res$targ_x), 3), labels = c('28', '33', '38')) #adding eccentricity
res$height <- factor(cut(res$targ_y, 3, labels = c('top', 'mid', 'bottom'))) #adding target height
#finally - target location
res$target <- paste0(res$ecc, res$height)
res$group <- factor(substr(res$subject_nr, 1, 1)) # adding group: 1 = HC, 2/4 = patient
res$site <- factor(substr(res$subject_nr, 1, 1)) # 1.2 edinburgh, 3/4 norwich
# group '4' = patients from norwich, change to group = 2
for (i in 1:length(res$group)){
  if (isTRUE(res$group[i] == "4")){
    res$group[i] = 2
  }
}
for (i in 1:length(res$site)){
  if (isTRUE(res$site[i] == "2")){
    res$site[i] = 1
  }
}
res$group <- factor(res$group)
res$site <- factor(res$site)

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
res <- res[which(res$eye_move == 0 & res$void == 0), c(1,11:15,2:5,18:23,6:8,16:17)]

# add demographic information to this data
patient_demos <- read.csv('patient_demographics.csv') #loading patient demographics
control_demos <- read.csv('control_demographics.csv') #loading control demos
#extracting ACE data into seperate data-frame
ACEscores <- patient_demos[ ,c(1, 8:13)]
#isolating patient demographic information to bind with control
patient_demos <- patient_demos[, c(1:6)]
demo <- rbind(control_demos, patient_demos)
#remove patient 409 - no lateral reaching data here
demo <- demo[demo$subject_nr != 409, ]

#merging demo with res medians
res <- merge(demo, res)

## renaming variables to match radial reaching
names(res)[1] <- 'PPT'
names(res)[6] <- 'DIAGNOSIS'
names(res)[7] <- 'VIEW'
names(res)[8] <- 'SIDE'
names(res)[9] <- 'POSITION'
names(res)[12] <- 'TARGx'
names(res)[13] <- 'TARGy'
names(res)[14] <- 'LANDx'
names(res)[15] <- 'LANDy'
names(res)[23] <- 'RT'
names(res)[24] <- 'MT'
names(res)[25] <- 'GRP'
names(res)[26] <- 'SITE'

setwd(anaPath)
write.csv(res, "lateral-reaching_compiled.csv", row.names = FALSE)

ggplot(res) + geom_point(aes(x = TARGx, y = TARGy), shape = 4, size = 3) +
  geom_point(aes(x = LANDx, y = LANDy, colour = POSITION), shape = 1, size = 2) +
  facet_wrap(. ~PPT*VIEW) 
ggsave('lateral-reach_err.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, width = 15, height = 10, units = 'in', path = anaPath)


######### data aggregation + plotting ############
res_medians <- aggregate(AEdeg ~ PPT * SIDE * VIEW * POSITION * SITE * GRP * DIAGNOSIS, median, data = res)
colnames(res_medians)[colnames(res_medians)=='AEdeg'] <- 'AEmed' #change name to be more logical

# changing levels to be more informative
res_medians$SIDE <- factor(res_medians$SIDE, levels = c('left', 'right'))
levels(res_medians$SIDE) <- c('Left', 'Right')
levels(res_medians$GRP) <- c('Control', 'Patient') #changing group name from 1 = control, 2 = AD
levels(res_medians$VIEW) <- c('Peripheral', 'Free')
levels(res_medians$SITE) <- c('UOE', 'UEA')
res_medians$DIAGNOSIS <- factor(res_medians$DIAGNOSIS)
res_medians <- res_medians[order(res_medians$PPT), ] 

res_means <- aggregate(AEmed ~ PPT * SIDE * VIEW * SITE * GRP * DIAGNOSIS, mean, data = res_medians)
res_means <- res_means[order(res_means$PPT), ] 
colnames(res_means)[colnames(res_means) == 'AEmed'] <- 'AEmean'
# save data
write.csv(res_medians, 'lateral-reaching_medians.csv', row.names = FALSE)
write.csv(res_means, 'lateral-reaching_means.csv', row.names = FALSE)

# to calculate PMI need to cast by task....
PMIdata <- dcast(res_means, PPT+GRP+SITE+SIDE+DIAGNOSIS ~ VIEW) #different data-frame
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
                    strip.text.x = element_text(size = 10)) -> eccPlot


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
                     strip.text.x = element_text(size = 10)) -> eccPlot

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
                     strip.text.x = element_text(size = 10)) -> eccPlot

ggsave('AD_ecc.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

##### PMI collapsed across side #####
PMIav <- aggregate(PMI ~ PPT * SITE * GRP * DIAGNOSIS, mean, data = PMIdata)
PMIav <- PMIav[order(PMIav$PPT), ]

##### outlier calculation, for controls only ######
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

###### directional error calc ######
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

########### next steps: comparing patients to controls
