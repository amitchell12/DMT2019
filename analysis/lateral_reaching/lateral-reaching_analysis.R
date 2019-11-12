library(readr)
library(ggplot2)
library(Rmisc)
library(gridExtra)
library(plyr)
library(tidyverse)

#set working directory to where data is
#on mac
#anaPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/analysis/lateral_reaching/control_data'
#dataPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/rawdata'
#on pc
dataPath <- 'S:/groups/DMT/data/control'
anaPath <- 'S:/groups/DMT/analysis/lateral_reaching/control_data'
setwd(dataPath)

#variable information
#nParticipants = 1:2 #for testing
nParticipants = 1:24 #for analysis
nTasks = 1:4 #total number of tasks - free + peripheral reaching, visual detection (non-dominant, dominant)
nSide = 1:2 #total number of sides tested - left, right
nTargets = 1:9
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

left_PMIdata = data_frame() #empty dataframe for saving all PP data
right_PMIdata = data_frame()

#make complete list of RESULTS files (trial information)
# for reaching only - visual detection analysed seperately
for (x in nParticipants){
  if (x < 10){ #getting around issue of added '0' on PP 1-9
    xx = sprintf('0%d', x)
  } else {
    xx = x
  }
  # creating paths
  setwd(dataPath)
  controlID = sprintf("1%s", xx)
  ppPath = (file.path(dataPath, controlID))
  setwd(anaPath)
  dir.create(controlID)
  
  resfiles = list.files(path=ppPath, full.names=TRUE, recursive = FALSE, 
                             include.dirs = FALSE, pattern = "*.csv")
  
  for (t in nTasks){
    setwd(dataPath) #make sure we are in this folder
    task_name = switch(t, "peripheral_left", "peripheral_right","free_left", "free_right")
    res <- read.csv(resfiles[t])
    # assigning key data points to dataframe for each task
    res1 = data.frame(res$targ_x, res$targ_y, res$land_x, res$land_y,
                      res$reach_duration, res$time_touch_offset, res$target_onset,
                      res$eye_move, res$void_trial)
    names(res1) <- gsub("res.", "", names(res1)) #renaming vectors without 'res.'
    
    # now to remove eye-move and void trials from each task
    res1 <- res1[res1$eye_move == 0, ]
    res1 <- res1[res1$void_trial == 0, ]

    # calculating x and y error for each targ location
    res1$xerr = res1$land_x - res1$targ_x
    res1$yerr = res1$land_y - res1$targ_y
    # transforming pixels into mm 
    res1$xerr_mm = res1$xerr*mm_perPix
    res1$yerr_mm = res1$yerr*mm_perPix
    # transforming pixels into degrees
    res1$xerr_deg = res1$xerr*deg_perpix
    res1$yerr_deg = res1$yerr*deg_perpix
    
    #absolute error in mm
    res1$abserr_mm = sqrt(res1$xerr_mm^2 + res1$yerr_mm^2)
    res1$abserr_deg = sqrt(res1$xerr_deg^2 + res1$yerr_deg^2)
    
    #group variables
    res1$ecc <- factor(cut(res1$targ_x, 3, label=c('28','33','38')))
    res1$height <- factor(cut(res1$targ_y, 3, label=c('top', 'mid', 'bottom')))
    res1$target <- paste0(res1$ecc, res1$height)
    
    file_name = sprintf("%s", task_name)
    assign(file_name, res1) #assigning csv file logical name, yay
    
    # saving filtered file as is to analysis folder for PP
    pp_anaPath = file.path(anaPath, controlID)
    setwd(pp_anaPath)
    write.csv(res1, sprintf("%s_allData.csv", file_name), row.names = TRUE)
    
    #median for each target location
    targ_meds <- res1 %>% #creating another data frame
      group_by(target) %>% 
      summarise(AEmm = median(abserr_mm), AEdeg = median(abserr_deg))
    
    #median for each eccentricity
    ecc_meds <- res1 %>%
      group_by(ecc) %>%
      summarise(AEmm = median(abserr_mm), AEdeg = median(abserr_deg), 
                AEmm_SD = sd(abserr_mm), AEdeg_SD = sd(abserr_deg))
    
    targ_meds$ecc <- substr(targ_meds$target, 1, 2) #adding column for eccentricity info only
    targ_meds$height <- substr(targ_meds$target, 3, 8) #adding column for height info only

    
    #adding name of task to data-frame for later compiling
    targ_meds$task <- file_name
    ecc_meds$task <- file_name
    #individually naming each data frame to task - to allow for compiling (median eccentricity)
    medfile_name = sprintf("%s_meds", task_name)
    assign(medfile_name, targ_meds)
    
    #individually naming each eccentricity median file 
    eccfile_name = sprintf("%s_eccmeds", task_name)
    assign(eccfile_name, ecc_meds)
    
  }
  
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
  
}

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
