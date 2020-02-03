##### COMPARING LATERAL AND RADIAL REACHING DATA #####
library(readr)
library(ggplot2)
library(reshape2)
library(ggpubr)

##### correlations #####
#anaPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/analysis' #mac
anaPath <- 'S:/groups/DMT/analysis' #pc
# loading lateral reaching data
#latPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/analysis/lateral_reaching'
latPath <- 'S:/groups/DMT/analysis/lateral_reaching'
setwd(latPath)
latData <- read.csv('lateral-reaching_PMI.csv')
# loading radial reaching data
#radPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/analysis/radial_reaching'
radPath <- 'S:/groups/DMT/analysis/radial_reaching'
setwd(radPath)
radData <- read.csv('radial-reaching_PMI.csv')
# TVA path
TVApath <- 'S:/groups/DMT/analysis/TVA/'
setwd(TVApath)
TVAData <- read.csv('TVA_values.csv')

# making sure all matrices are the same
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
dat_side <- rbind(latData, radData)
corrDat_side <- dcast(dat_side, PPT+GRP+SIDE ~ TASK, value.var = 'PMI')
dat <- aggregate(PMI ~ TASK * GRP * PPT, mean, data = dat_side) 
corrDat <- dcast(dat, PPT+GRP ~ TASK, value.var = 'PMI')

ggscatter(corrDat_side, x = 'Lateral', y = 'Radial', add = 'reg.line', conf.int = TRUE,
          cor.coef = TRUE, cor.method = 'pearson') + 
  facet_grid(cols = vars(SIDE), rows = vars(GRP)) + 
  ylab('Radial reaching PMI (deg)') + xlab('Lateral reaching PMI (deg)')

ggsave('task_correlations_pearson.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

#file with data organised by participant
names(latData)[4] <- 'Free_lat'
names(latData)[5] <- 'Perip_lat'
names(latData)[6] <- 'PMI_lat'
names(radData)[4] <- 'Free_rad'
names(radData)[5] <- 'Perip_rad'
names(radData)[6] <- 'PMI_rad'
allDat <- merge(latData, radData, by = c('PPT', 'GRP', 'SIDE'), all = TRUE)
allDat <- allDat[c(1:5, 8:9, 6, 10) ]

###### TVA data ######
# ranaming to match reaching data-frame
names(TVAData)[2] <- 'PPT'
TVAData <- TVAData[c(2:7)]
# merge with allDat
allDat <- merge(allDat, TVAData, by = c('PPT', 'GRP'), all = TRUE)

setwd(anaPath)
write.csv(allDat, 'all-task_results.csv', row.names = FALSE)

