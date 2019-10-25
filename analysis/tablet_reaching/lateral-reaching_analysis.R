library(readr)
library(ggplot2)
library(Rmisc)
library(gridExtra)
library(plyr)

#set working directory to where data is
#on mac
anaPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/analysis/lateral_reaching'
dataPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/rawdata'
#on pc
#dataPath <- 'S:/groups/DMT/data/control'
#anaPath <- 'S:/groups/DMT/analysis/lateral_reaching'
setwd(dataPath)

#variable information
nParticipants = 1 #for testing
#nParticipants = 1:24 #for analysis
nTasks = 1:4 #total number of tasks - free + peripheral reaching, visual detection (non-dominant, dominant)
nSide = 1:2 #total number of sides tested - left, right
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
    
  
    
    file_name = sprintf("%s", task_name)
    assign(file_name, res1) #assigning csv file logical name, yay
    
    # saving filtered file as is to analysis folder for PP
    pp_anaPath = file.path(anaPath, ppID)
    setwd(pp_anaPath)
    write.csv(res1, sprintf("%s.csv", file_name), row.names = TRUE)
  }
  
  #calculating absolute error for each target location
  
  
}
  

#make empty dataframe to stack up trial information
#CLF <- read.csv(text = "avg_rt,count_committed,count_valid_trial,datetime,eye_move,fix_x,land_x,land_y,reach_duration,response_time_target_mouse_response,sanity_check,subject_nr,subject_parity,targ_x,targ_y,target_onset,targets_file,time_touch_offset,touch_x")


#stack them all up
#for(x in resfiles){
#  tmp <- read.csv(x)
#}
  

CLF$x_Err <- CLF$land_x-CLF$targ_x
CLF$y_Err <- CLF$land_y-CLF$targ_y

CLF[is.na(CLF$eye_move), "eye_move"] <- 0

CLF$ECC <- factor(cut(CLF$targ_x, 3, label=c("N", "M", "F")))
CLF <- CLF[CLF$eye_move == 0, ]

CLF$HEIGHT <- factor(cut(CLF$targ_y, 3, label=c("L", "M", "H")))

CLF$TARGET <- paste0(CLF$ECC, CLF$HEIGHT)

CLF$ABS_ERR <- sqrt(CLF$x_Err^2 + CLF$y_Err^2)


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
