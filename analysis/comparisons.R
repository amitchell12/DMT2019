##### COMPARING LATERAL AND RADIAL REACHING DATA #####
library(readr)
library(ggplot2)
library(reshape2)
library(ggpubr)

##### correlations #####
#anaPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/analysis' #mac
anaPath <- 'S:/groups/DMT/analysis/task-correlations' #pc
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
          cor.coef = TRUE, cor.method = 'spearman') + 
  facet_grid(cols = vars(SIDE), rows = vars(GRP)) + 
  ylab('Radial reaching PMI (deg)') + xlab('Lateral reaching PMI (deg)')

ggsave('reaching_correlations_spearman.png', plot = last_plot(), device = NULL, dpi = 300, 
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

# removing NA values for plotting
corr_allDat <- na.omit(allDat)

##### TVA lateral reaching
### correlate K and lateral reaching
ggscatter(corr_allDat, x = 'PMI_lat', y = 'K', add = 'reg.line', conf.int = TRUE,
          cor.coef = TRUE, cor.method = 'spearman') + 
  facet_grid(cols = vars(SIDE), rows = vars(GRP)) + 
  ylab('vSTM') + xlab('Lateral reaching PMI (deg)')

ggsave('LAT-k_side.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

# collapsed across side
KLAT <- aggregate(PMI_lat ~ K * PPT * GRP, mean, data = corr_allDat)
ggscatter(KLAT, x = 'PMI_lat', y = 'K', add = 'reg.line', conf.int = TRUE,
          cor.coef = TRUE, cor.method = 'spearman') + 
  facet_wrap(~GRP) + 
  ylab('vSTM') + xlab('Lateral reaching PMI (deg)')

ggsave('LAT-k_grp.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

### correlate C and lateral reaching
ggscatter(corr_allDat, x = 'PMI_lat', y = 'C', add = 'reg.line', conf.int = TRUE,
          cor.coef = TRUE, cor.method = 'spearman') + 
  facet_grid(cols = vars(SIDE), rows = vars(GRP)) + 
  ylab('Processing speed (ms)') + xlab('Lateral reaching PMI (deg)')

ggsave('LAT-C_side.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

# collapsed across side
CLAT <- aggregate(PMI_lat ~ C * PPT * GRP, mean, data = corr_allDat)
ggscatter(CLAT, x = 'PMI_lat', y = 'C', add = 'reg.line', conf.int = TRUE,
          cor.coef = TRUE, cor.method = 'spearman') + 
  facet_wrap(~GRP) + 
  ylab('Processing speed (ms)') + xlab('Lateral reaching PMI (deg)')

ggsave('LAT-C_grp.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

### correlate K and radial reaching
ggscatter(corr_allDat, x = 'PMI_rad', y = 'K', add = 'reg.line', conf.int = TRUE,
          cor.coef = TRUE, cor.method = 'spearman') + 
  facet_grid(cols = vars(SIDE), rows = vars(GRP)) + 
  ylab('vSTM') + xlab('Radial reaching PMI (deg)')

ggsave('RAD-k_side.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

# collapsed across side
KRAD <- aggregate(PMI_rad ~ K * PPT * GRP, mean, data = corr_allDat)
ggscatter(KRAD, x = 'PMI_rad', y = 'K', add = 'reg.line', conf.int = TRUE,
          cor.coef = TRUE, cor.method = 'spearman') + 
  facet_wrap(~GRP) + 
  ylab('vSTM') + xlab('Lateral reaching PMI (deg)')

ggsave('RAD-k_grp.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

### correlate C and radial reaching
ggscatter(corr_allDat, x = 'PMI_rad', y = 'C', add = 'reg.line', conf.int = TRUE,
          cor.coef = TRUE, cor.method = 'spearman') + 
  facet_grid(cols = vars(SIDE), rows = vars(GRP)) + 
  ylab('Processing speed (ms)') + xlab('Radial reaching PMI (deg)')

ggsave('RAD-C_side.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

# collapsed across side
CRAD <- aggregate(PMI_rad ~ C * PPT * GRP, mean, data = corr_allDat)
ggscatter(CRAD, x = 'PMI_rad', y = 'C', add = 'reg.line', conf.int = TRUE,
          cor.coef = TRUE, cor.method = 'spearman') + 
  facet_wrap(~GRP) + 
  ylab('Processing speed (ms)') + xlab('Lateral reaching PMI (deg)')

ggsave('RAD-C_grp.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)


