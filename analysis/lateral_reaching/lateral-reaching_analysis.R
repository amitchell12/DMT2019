library(readr)
library(ggplot2)
library(Rmisc)
library(gridExtra)
library(plyr)
library(tidyverse)

#set working directory to where data is
#on mac
#anaPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/analysis/lateral_reaching'
#dataPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/rawdata'
#on pc
dataPath <- 'S:/groups/DMT/data/control'
anaPath <- 'S:/groups/DMT/analysis/lateral_reaching'
setwd(dataPath)

#variable information
nParticipants = 2 #for testing
#nParticipants = 1:24 #for analysis
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

#make complete list of RESULTS files (trial information)
# for reaching only - visual detection analysed seperately
for (x in nParticipants){
  # creating paths
  ppID = sprintf("10%s", x)
  ppPath = (file.path(dataPath, ppID))
  setwd(anaPath)
  dir.create(ppID)
  
  resfiles = list.files(path=ppPath, full.names=TRUE, recursive = FALSE, 
                             include.dirs = FALSE, pattern = "*.csv")
  
  for (t in nTasks){
    task_name = switch(t, "peripheral_left", "peripheral_right","free_left", "free_right",
                       "detection_left", "detection_right")
    res <- read.csv(resfiles[t])
    # assigning key data points to dataframe for each task
    res1 = data.frame(res$targ_x, res$targ_y, res$land_x, res$land_y,
                      res$reach_duration, res$time_touch_offset, res$target_onset,
                      res$eye_move, res$void_trial)
    names(res1) <- gsub("res.", "", names(res1)) #renaming vectors without 'res.'
    
    # now to remove eye-move and void trials from each task
    res1 <- res1[res1$eye_move == 0, ]
    res1 <- res1[res1$void_trial == 0, ]
    
    #adding extra info columns
    
    
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
    pp_anaPath = file.path(anaPath, ppID)
    setwd(pp_anaPath)
    write.csv(res1, sprintf("%s_allData.csv", file_name), row.names = TRUE)
    
    #median for each target location
    targ_meds <- res1 %>% #creating another data frame
      group_by(target) %>% 
      summarise(meds_mm = median(abserr_mm), meds_deg = median(abserr_deg))
    
    targ_meds$ecc <- substr(targ_meds$target, 1, 2) #adding column for eccentricity info only
    targ_meds$height <- substr(targ_meds$target, 3, 8) #adding column for height info only

    
    #adding name of task to data-frame for later compiling
    targ_meds$task <- file_name
    #individually naming each data frame to task - to allow for compiling
    medfile_name = sprintf("%s_meds", task_name)
    assign(medfile_name, targ_meds) 
    
    #plot the medians - then save plot to each PP folder
    ggplot(targ_meds, aes(x = ecc, y = meds_deg, colour = target)) +
      geom_point(size = 4) + ylim(0,4) + 
      geom_hline(yintercept = 1, linetype = 'dotted') +
      theme_bw()
    
  }

  
  #append each median task to each other then calculate means for each task
  
  
  
}
  


#boxplots
ggplot(CLF, aes(x=SUB, y=ABS_ERR, colour=ECC)) +
  geom_boxplot(outlier.alpha=0) +
  geom_jitter(position=position_dodge(width=.8), size=4, alpha=.5) +
  ylim(0,300) + geom_hline(yintercept=55, linetype="dotted")+
  theme_bw()


#mean of medians
mean_AE <- summarySE(data=CLF, measurevar = "ABS_ERR", groupvars = c("SUB", "GRP"))

ggplot(mean_AE, aes(x = SUB, y=ABS_ERR, colour=GRP)) + geom_point(size=5, alpha=.5) +
  geom_errorbar(aes(ymin=ABS_ERR-ci, ymax=ABS_ERR+ci)) + ylim(c(50,175))-> ALL

grid.arrange(ALL, NOTNEAR, ncol=2)

mean_CE <- summarySE(data=CLF, measurevar = "x_Err", groupvars = c("SUB", "GRP"))

ggplot(mean_CE, aes(x = SUB, y=x_Err, colour=GRP)) + geom_point(size=5, alpha=.5) +
  geom_errorbar(aes(ymin=x_Err-ci, ymax=x_Err+ci))

lm_fits <- read.csv(text="SUB,rsq,a,b")

row=1


for(SUB in levels(CLF$SUB)){
  tmp <- CLF[CLF$SUB==SUB, ]
  model <- lm(ABS_ERR~ECC, data=tmp)
  lm_fits[row, 1] <- SUB
  lm_fits[row, 3:4] <- coefficients(model)[1:2]
  lm_fits[row, 2] <- summary(model)$r.squared
  row=row+1
}


med <- aggregate(ABS_ERR~SUB*ECC, median, data=CLF[CLF$targ_x > 300, ])

mean_med <- aggregate(ABS_ERR~SUB, mean, data=med)

ggplot(mean_med, aes(x=SUB, y=ABS_ERR))+ geom_point(size=5)
