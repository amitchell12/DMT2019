##### COMPARING LATERAL AND RADIAL REACHING DATA #####
library(readr)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(plyr)

##### correlations #####
## PATHS WILL NEED TO CHANGE IF USING DIFF COMPUTER ##
anaPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/analysis/task-correlations' #mac
# loading lateral reaching data
latPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/analysis/lateral_reaching'
setwd(latPath)
latData <- read.csv('lateralPMI_all.csv')
# loading radial reaching data
radPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/analysis/radial_reaching'
setwd(radPath)
radData <- read.csv('radialPMI.csv')

# adding extra info necessary to merge
latData$TASK <- 'Lateral'
radData$TASK <- 'Radial'
# changing labels of 'left/right' in rad data so they match lateral
latData$SIDE <- factor(latData$SIDE, labels = c('Left', 'Right'))
latData <- latData[, c(1:7,9:12)]
# labelling groups in rad data
radData$GRP <- factor(radData$GRP)
# revalue UEA numbering system to match EDB
radData$GRP <- revalue(radData$GRP, c("3"="1", "4"="2"))
radData$GRP <- factor(radData$GRP, labels = c('Control','Patient'))

##### PMI CORRELATIONS ######
allDat <- rbind(latData, radData)
# average across side
allDat <- aggregate(PMI ~ TASK*PPT*GRP*DIAGNOSIS*SITE*AGE, mean, data = allDat)
# cast across task
corrPMI <- dcast(allDat, PPT+GRP+DIAGNOSIS ~ TASK, value.var = 'PMI')
corrPMI$DIAGNOSIS <- factor(corrPMI$DIAGNOSIS, levels = c('HC','MCI','AD'))

## running correlation - split between patients
corrPMI_patients <- corrPMI[corrPMI$GRP == 'Patient' ,]
corrPMI_controls <- corrPMI[corrPMI$GRP == 'Control' ,]

res_control <- cor.test(corrPMI_controls$Lateral, corrPMI_controls$Radial, 
                method = "spearman")
res_patients <- cor.test(corrPMI_patients$Lateral, corrPMI_patients$Radial, 
                        method = "spearman")

ggplot(corrPMI, aes(x = Lateral, y = Radial)) + 
  geom_point(aes(shape = DIAGNOSIS, colour = DIAGNOSIS), size = 2.5) +
  geom_smooth(method=lm, se=FALSE, colour = 'black') +
  scale_colour_manual(values = c('black','grey30','grey60')) +
  facet_wrap(~GRP) + 
  xlim(0,50) +
  ylab('Radial reaching PMI (mm)') + xlab('Lateral reaching PMI (mm)') +
  theme_classic() + theme(legend.position = 'bottom',
                          legend.title = element_blank(),
                          legend.text = element_text(size = 10),
                          axis.text = element_text(size = 10),
                          axis.title = element_text(size = 12),
                          strip.text = element_text(size = 12))

ggsave('reach_correlations.png', plot = last_plot(), device = NULL, dpi = 300, 
       width = 8, height = 5, path = anaPath)
