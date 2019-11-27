##### comparing young to old
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
dataPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/analysis/pointing'
anaPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/analysis'
setwd(dataPath)

########### variable info ###########
visAngle <- function(size, distance){
  # this function calculates visual angle
  # size and distance must be in the same units
  Rad = 2*atan(size/(2*distance))
  Ang = Rad*(180/pi)
  return(Ang)
}

############ file info ############
#getting all datafiles and compiling (patient + control)
filenames <- dir(dataPath, recursive = TRUE, full.names = FALSE, pattern = '.csv')

# making results file
res <- read.csv(text = "subject_nr,targ_x,targ_y,land_x,land_y,task")

for (file in filenames){
  if (isTRUE(substr(basename(file), 8, 8)=='9')){
  tmp <- read.csv(file)[, c(12,14:15,7:8)]
  tmp$task <- 'periph'
  
  res<- rbind(tmp, res)
  }
}

# adding key details to data-frame
res$task <- factor(res$task) # task = factor
res$task = revalue(res$task, c('i'='periph', 'r'='free')) #renaming task
res$side <- factor(res$targ_x < 0, label = c('right', 'left')) #adding factor of side
res$ecc <- factor(cut(abs(res$targ_x), 3), labels = c('28', '33', '38')) #adding eccentricity
res$height <- factor(cut(res$targ_y, 3, labels = c('top', 'mid', 'bottom'))) #adding target height
#finally - target location
res$target <- paste0(res$ecc, res$height)

# calculating x and y error for each targ location
res$xerr_mm = (res$land_x - res$targ_x)*mm_perPix # in mm
res$yerr_mm = (res$land_y - res$targ_y)*mm_perPix
res$xerr_deg = visAngle(size= res$xerr_mm, distance= 400) # in deg
res$yerr_deg = visAngle(size= res$yerr_mm, distance= 400)

#absolute error in mm
res$AEmm = sqrt(res$xerr_mm^2 + res$yerr_mm^2)
res$AEdeg = sqrt(res$xerr_deg^2 + res$yerr_deg^2)

######### data aggregation + plotting ############
res_medians <- aggregate(AEdeg ~ ecc * side * task * subject_nr, median, data = res)
colnames(res_medians)[colnames(res_medians)=='AEdeg'] <- 'AEmed' #change name to be more logical
res_means <- aggregate(AEmed ~ task * side * subject_nr, mean, data = res_medians)
colnames(res_means)[colnames(res_means) == 'AEmed'] <- 'AEmean'

## load lateral reaching means to plot and compare
latPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/analysis/lateral_reaching'
setwd(latPath)
# load data
lat_means <- read.csv('lateral-reaching_means.csv')
lat_means <- lat_means[lat_means$group == '1', c(1:3,5)]
lat_means <- lat_means[lat_means$task == 'periph', ]
colnames(lat_means)[colnames(lat_means) == 'AEdeg'] <- 'AEmean'

#compiling data
all_means <- rbind(lat_means, res_means)
all_means$group <- factor(substr(all_means$subject_nr, 1, 1))
levels(all_means$group) <- c('Elderly', 'Young')

# plotting periphral reaching data
ggplot(all_means, aes(x = side, y = AEmean, colour = group)) +
  geom_point(shape = 1, size = 2) +
  geom_line(aes(group = subject_nr), size = 0.5, alpha = .5) +
  facet_wrap(~group) + ylim(-.5,5) +
  labs(x = 'Side', y = 'Mean AE (deg)', element_text(size = 12)) +
  scale_colour_manual(values = c('black', 'grey50')) +
  stat_summary(aes(y = AEmean, group = 1), fun.y = mean, colour = "red", 
               geom = 'point', shape = 3, stroke = 1, size = 2, group = 1, alpha = .7) +
  theme_bw() + theme(legend.position = 'none', text = element_text(size = 10),
                     strip.text.x = element_text(size = 8)) -> meansPlot

ggsave('young-old_plot.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)



