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

# calculating x and y error for each targ location
res$xerr_mm = (res$land_x - res$targ_x)*mm_perPix # in mm
res$yerr_mm = (res$land_y - res$targ_y)*mm_perPix
res$xerr_deg = visAngle(size= res$xerr_mm, distance= 400) # in deg
res$yerr_deg = visAngle(size= res$yerr_mm, distance= 400)

#absolute error in mm
res$AE = sqrt(res$xerr_mm^2 + res$yerr_mm^2)
res$AEdeg = sqrt(res$xerr_deg^2 + res$yerr_deg^2)

# add demographic information to this data
patient_demos <- read.csv('patient_demographics.csv') #loading patient demographics
control_demos <- read.csv('control_demographics 2.csv') #loading control demos
#extracting ACE data into seperate data-frame
ACEscores <- patient_demos[ ,c(1, 8:13)]
#isolating patient demographic information to bind with control
patient_demos <- patient_demos[, c(1:6)]
demo <- rbind(control_demos, patient_demos)

#merging demo with res medians
res <- merge(demo, res, by = 'subject_nr')

# eye move and void trials
# counting eye-move per participant
setwd(anaPath)
nEye_move <- aggregate(eye_move ~ subject_nr * diagnosis * task, sum, data = res)
nTrials <- res %>% 
  group_by(subject_nr,task) %>% 
  tally()

# getting percentage eye-move for each group, and each condition
Eye_move <- merge(nEye_move, nTrials)
sumEye_move <- aggregate(eye_move ~ diagnosis * task, sum, data = Eye_move)
sumTrials <- aggregate(n ~ diagnosis * task, sum, data = Eye_move)
totEye_move <- merge(sumEye_move, sumTrials)

totEye_move$per <- (totEye_move$eye_move/totEye_move$n)*100

write.csv(totEye_move, 'eye-movements.csv', row.names = FALSE)

nVoid <- aggregate(void_trial ~ subject_nr * diagnosis * task, sum, data = res)
tot_void <- aggregate(void_trial ~ diagnosis * task, sum, data = nVoid)

# removing and reorganising
res <- res[which(res$eye_move == 0 & res$void == 0), c(1:6,16:19,7:13,23:28,21,22)]

## renaming variables to match radial reaching
names(res)[1] <- 'PPT'
names(res)[3] <- 'AGE'
names(res)[4] <- 'HAND'
names(res)[5] <- 'ED'
names(res)[6] <- 'DIAGNOSIS'
names(res)[7] <- 'VIEW'
names(res)[8] <- 'SIDE'
names(res)[9] <- 'POSITION'
names(res)[11] <- 'TARGx'
names(res)[12] <- 'TARGy'
names(res)[13] <- 'LANDx'
names(res)[14] <- 'LANDy'
names(res)[16] <- 'RT'
names(res)[17] <- 'MT'
names(res)[24] <- 'GRP'
names(res)[25] <- 'SITE'

## getting dominant + non-dominant sides, for analysis
res$HAND <- factor(res$HAND, labels = c('right','left'))
# if hand = side, dominant; else non dominant
res$DOM <- as.numeric(res$SIDE == res$HAND) #1 = dominant, 0 = non-dominant
res$DOM <- factor(res$DOM, labels= c('ND','D'))
# change order so dominance up-front
res <- res[, c(1:8,26,9:25)]

## remove trial with wildly high RT (> 6000ms)
res <- res[res$RT < 6000 ,]

write.csv(res, "lateral-reaching_compiled.csv", row.names = FALSE)

