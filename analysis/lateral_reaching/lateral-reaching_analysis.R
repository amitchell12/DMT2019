library(readr)
library(ggplot2)
library(Rmisc)
library(gridExtra)
library(plyr)
library(tidyverse)
library(reshape2)

#set working directory to where data is
#on mac
anaPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/analysis/lateral_reaching'
dataPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/rawdata'
#on pc
#dataPath <- 'S:/groups/DMT/data'
#anaPath <- 'S:/groups/DMT/analysis/lateral_reaching'
setwd(dataPath)

########### variable info ###########
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

# removing df (subject number = 300)
res <- res[res$subject_nr < 300, ]

# adding key details to data-frame
res$task <- factor(res$task) # task = factor
res$task = revalue(res$task, c('i'='periph', 'r'='free')) #renaming task
res$side <- factor(res$targ_x < 0, label = c('right', 'left')) #adding factor of side
res$ecc <- factor(cut(abs(res$targ_x), 3), labels = c('28', '33', '38')) #adding eccentricity
res$height <- factor(cut(res$targ_y, 3, labels = c('top', 'mid', 'bottom'))) #adding target height
#finally - target location
res$target <- paste0(res$ecc, res$height)
res$group <- factor(substr(res$subject_nr, 1, 1)) # adding group: 1 = HC, 2 = patient

# eye move and void trials
# counting eye-move per participant
nEye_move <- aggregate(res$eye_move, by=list(subject_nr = res$subject_nr), FUN=sum)
nVoid <- aggregate(res$void, by=list(subject_nr = res$subject_nr), FUN=sum)

# calculating x and y error for each targ location
res$xerr_mm = (res$land_x - res$targ_x)*mm_perPix # in mm
res$yerr_mm = (res$land_y - res$targ_y)*mm_perPix
res$xerr_deg = (res$land_x - res$targ_x)*deg_perpix # in deg
res$yerr_deg = (res$land_y - res$targ_y)*deg_perpix

#absolute error in mm
res$AEmm = sqrt(res$xerr_mm^2 + res$yerr_mm^2)
res$AEdeg = sqrt(res$xerr_deg^2 + res$yerr_deg^2)

# removing and reorganising
res <- res[which(res$eye_move == 0 & res$void == 0), c(1,11,2:5,17:22,6:8,12:16)]

setwd(anaPath)
write.csv(res, "lateral-reaching_compiled.csv", row.names = FALSE)

ggplot(res) + geom_point(aes(x = targ_x, y = targ_y), shape = 4, size = 3) +
  geom_point(aes(x = land_x, y = land_y, colour = ecc), shape = 1, size = 2) +
  facet_wrap(. ~subject_nr*task) -> allPP_plot

######### data aggregation + plotting ############
res_medians <- aggregate(AEdeg ~ ecc * side * subject_nr * group, median, data = res)
res_means <- aggredate(AEdeg ~ ecc * side * subject_nr * group)

ggplot(res_medians) + geom_point(aes(x = ecc, y = AEdeg, colour = subject_nr)) +
  facet_wrap(~side)

left_PMIdata = data_frame() #empty dataframe for saving all PP data
right_PMIdata = data_frame()



  
  setwd(pp_anaPath) #make sure we are in this folder
  #append free-peripheral task-med frames for each side
  # median for each eccentricity
  left_data <- rbind(free_left_eccmeds, peripheral_left_eccmeds)
  left_data$side <- 'Left'
  right_data <- rbind(free_right_eccmeds, peripheral_right_eccmeds)
  right_data$side <- "Right"
  median_data <- rbind(left_data, right_data)
  # saving these
  write.csv(median_data, sprintf("%s_AEmedians.csv", controlID), row.names = TRUE) 

  median_data$condition <- substr(median_data$task, 1, 4)

  
  ## calculating overall mean for each side
  meanAE_mm = summarySE(data = median_data, measurevar = "AEmm", groupvars = c("task", "side"))
  meanAE_deg = summarySE(data = median_data, measurevar = "AEdeg", groupvars = c("task", "side"))
  
  meanAE = meanAE_deg
  # pp information to left mean data
  meanAE$sub <- controlID
  # saving
  write.csv(meanAE, sprintf("%s_meanAE.csv", controlID), row.names = TRUE) 
  
  # plotting! (each participant)
  # mean for each eccentricity
  plot_name = sprintf('%s_medians.png', controlID)
  means <- ggplot(median_data, aes(x = ecc, y = AEdeg, colour = condition)) +
    geom_point(size = 3, position = position_dodge(.5)) + ylim(-1,8) + 
    geom_errorbar(aes(ymin = AEdeg-AEdeg_SD, ymax = AEdeg+AEdeg_SD), width = .2, 
                  size = 0.7, position = position_dodge(.5)) +
    geom_hline(yintercept = 1, linetype = 'dotted') + 
    scale_colour_manual(values = c('free'='grey60', 'peri'='black'), labels = c('free', 'peripheral')) +
    labs(x = 'Target eccentricity (deg)', y = 'Absolute error (deg)') + 
    theme_bw() + theme(legend.position = 'bottom', legend.title = element_blank(),
                       legend.box = 'vertical')
  
  means + facet_grid(cols = vars(side)) + theme(strip.text.x = element_text(size = 11))
  plotPath = file.path(anaPath, 'plots')
  ggsave(plot_name, plot = last_plot(), device = NULL, dpi = 300, width = 7, height = 4,
         path = plotPath)
  
  # calculating PMI score - then plotting :)
  left = meanAE$AEdeg[3]-meanAE$AEdeg[1]
  right = meanAE$AEdeg[4]-meanAE$AEdeg[2]
  PMI <- data_frame(c(left, right))
  PMI$sub <- controlID
  PMI$side <- c('Left', 'Right')
  colnames(PMI)[1] <- 'PMI'

  #PMIname = sprintf("%s_PMI", controlID)
  #assign(PMIname, PMI)
  
  left_PMIdata <- rbind(left_PMIdata, PMI[1,]) #creating big data-frames
  right_PMIdata <- rbind(right_PMIdata, PMI[2,])
  


#compile all participant data
all_PMI = rbind(left_PMIdata, right_PMIdata)
all_PMI$side <- factor(all_PMI$side)
#summary data
all_PMIsummary <- summarySE(all_PMI, measurevar = "PMI", groupvars = "side")

# plot left and right :)
# need to get error bars on to this somehow - leave for now
plot_name = 'PMI_allPP.png'
PMIplot <- ggplot(all_PMI, aes(x = side, y = PMI, colour = sub)) +
  geom_point(size = 3, position = position_dodge(.1)) + ylim(-1,8) + 
  stat_summary(aes(y = PMI, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', group = 1) + 
  labs(x = '', y = 'Peripheral misreaching index (deg)', element_text(size = 13)) +
  theme_bw() + theme(legend.position = "right", legend.title = element_blank())

ggsave(plot_name, plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, width = 5, height = 6.5, path = plotPath)

plot_name = 'PMI_boxplot.png'
PMIbox <- ggplot(all_PMI, aes(x = side, y = PMI)) +
  geom_boxplot(size = 0.75, outlier.color = "black", outlier.shape = 16, 
               outlier.size = 2, notch = FALSE) +
  stat_summary(fun.y = mean, geom = 'point', shape = 18, size = 5) +
  ylim(-1,6) + 
  labs(x = '', y = 'Peripheral misreaching index (deg)', element_text(size = 13)) +
  theme_classic() + theme(legend.position = "right", legend.title = element_blank())

ggsave(plot_name, plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, width = 5, height = 6.5, path = plotPath)

## outlier calculation
# median values for each side
tmp <- aggregate(all_PMI$PMI,  by=list(side = all_PMI$side), FUN=median)
names(tmp)[2] <- 'med'
all_PMI <- merge(tmp, all_PMI)

# calculating MAD for each (absolute value)
all_PMI$AD <- abs(all_PMI$PMI - all_PMI$med)
tmp <- aggregate(all_PMI$AD,  by=list(side = all_PMI$side), FUN=median)
names(tmp)[2] <- 'MAD'
all_PMI <- merge(all_PMI, tmp)

# adjusted z-score from these values
all_PMI$az <- (all_PMI$PMI - all_PMI$med)/(all_PMI$MAD * 1.4826)
all_PMI$z <- scale(all_PMI$PMI)

plot_name = 'Zscore.png'
Zplot <- ggplot(all_PMI, aes(x = side, y = z, colour = sub)) +
  geom_point(size = 3, position = position_dodge(.1)) + ylim(-3,6) + 
  stat_summary(aes(y = z, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', group = 1) + 
  labs(x = '', y = 'Z-score', element_text(size = 13)) +
  theme_bw() + theme(legend.position = "right", legend.title = element_blank())

ggsave(plot_name, plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, width = 5, height = 6.5, path = plotPath)

plot_name = 'adjustedZ.png'
AZplot <- ggplot(all_PMI, aes(x = side, y = az, colour = sub)) +
  geom_point(size = 3, position = position_dodge(.1)) + ylim(-3,6) + 
  stat_summary(aes(y = az, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', group = 1) + 
  labs(x = '', y = 'Adjusted z-score', element_text(size = 13)) +
  theme_bw() + theme(legend.position = "right", legend.title = element_blank())

ggsave(plot_name, plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, width = 5, height = 6.5, path = plotPath)


sidecorr <- ggplot(all_PMI(x = ))