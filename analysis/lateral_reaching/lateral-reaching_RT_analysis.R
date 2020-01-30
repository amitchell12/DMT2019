##### analysing reaction time in young and old participants for the lateral reaching task
#### A.G.Mitchell 30.01.20
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
dataPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/analysis/lateral_reaching'
anaPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/analysis/lateral_reaching/lifespan/'
#on pc
#dataPath <- 'M:/Alex_Files/Experiments/DMT2019/raw_data/healthy'
#anaPath <- 'S:/groups/DMT/analysis/lateral-reaching/lifespan/'
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
mm_perPix = 1/pixPer_mm;############ file info ############

#getting file with key control and patient data
datfile <- read.csv('lateral-reaching_compiled.csv')
#ranaming
datfile$group <- factor(datfile$group)
levels(datfile$group) <- c('OA', 'AD')
levels(datfile$side) <- c('Left', 'Right')

##### looking reaction time data #####
# aggregating by trial
MT <- aggregate(reach_duration~group * ecc * task * side * subject_nr, mean, 
                                data= datfile)
RT <- aggregate(time_touch_offset~group * ecc * task * side * subject_nr, mean, 
               data= datfile)
ERR <- aggregate(AEdeg~group * ecc * task * side * subject_nr, median, 
                data= datfile)
# save all these
setwd(anaPath)
write.csv(MT, 'movement-time.csv', row.names = FALSE)
write.csv(RT, 'touch-offset-time.csv', row.names = FALSE)
write.csv(ERR, 'accuracy.csv', row.names = FALSE)


