##### COMPARING LATERAL AND RADIAL REACHING DATA #####
library(readr)
library(ggplot2)
library(reshape2)
library(ggpubr)

##### correlations #####
anaPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/analysis'
# loading lateral reaching data
latPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/analysis/lateral_reaching'
setwd(latPath)
latData <- read.csv('lateral-reaching_PMI.csv')
# loading radial reaching data
radPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/analysis/radial_reaching'
setwd(radPath)
radData <- read.csv('radial-reaching_PMI.csv')

# making sure both matrices are the same
names(latData)[1] <- 'PPT'
names(latData)[2] <- 'GRP'
names(latData)[3] <- 'SIDE'
names(latData)[4] <- 'Peripheral'
names(latData)[5] <- 'Free'
latData <- latData[c(1,2,3,5,4,6)] #reorgansing to fit radData

# adding extra info necessary to merge
latData$TASK <- 'Lateral'
radData$TASK <- 'Radial'
# merge
dat <- rbind(latData, radData)
corrdat <- dcast(dat, PPT+GRP+SIDE ~ TASK, value.var = 'PMI')

ggscatter(corrdat, x = 'Lateral', y = 'Radial', add = 'reg.line', conf.int = TRUE,
          cor.coef = TRUE, cor.method = 'pearson') + 
  facet_wrap(~GRP) + ylab('Radial reaching PMI (deg)') + xlab('Lateral reaching PMI (deg)')

ggsave('task_correlations.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)



