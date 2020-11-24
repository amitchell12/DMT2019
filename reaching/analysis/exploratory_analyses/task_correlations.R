##### COMPARING LATERAL AND RADIAL REACHING DATA #####
library(readr)
library(ggplot2)
library(reshape2)
library(ggpubr)

##### correlations #####
#anaPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/analysis/task-correlations' #mac
#anaPath <- 'S:/groups/DMT/analysis/task-correlations' #pc
anaPath <- "/Users/Alex/Documents/DMT/analysis/task-correlations/"
# loading lateral reaching data
#latPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/analysis/lateral_reaching'
#latPath <- 'S:/groups/DMT/analysis/lateral_reaching'
latPath <- "/Users/Alex/Documents/DMT/analysis/lateral_reaching/"
setwd(latPath)
latData <- read.csv('lateralPMI-filtered.csv')
# loading radial reaching data
#radPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/analysis/radial_reaching'
#radPath <- 'S:/groups/DMT/analysis/radial_reaching'
radPath <- "/Users/Alex/Documents/DMT/analysis/radial_reaching/"
setwd(radPath)
radData <- read.csv('radialPMI-filtered.csv')
# TVA path
#TVApath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/analysis/TVA/all/'
TVApath <- "/Users/Alex/Documents/DMT/analysis/TVA/all/"
setwd(TVApath)
TVAData <- read.csv('tva_values.csv')

# adding extra info necessary to merge
latData$TASK <- 'Lateral'
radData$TASK <- 'Radial'
latData_dir$TASK <- 'Lateral'
radData_dir$TASK <- 'Radial'
# changing labels of 'left/right' in rad data so they match lateral
latData$SIDE <- factor(latData$SIDE, labels = c('Left', 'Right'))
latData_dir$SIDE <- factor(latData_dir$SIDE, labels = c('Left', 'Right'))

## not now
#adding directional PMI into data-frame (along x-axis only)
#latData$dPMI <- latData_dir$PMI
#radData$dPMI <- radData_dir$PMI

##### PMI data ######
latData <- latData[, c(1:7,9:12)]

dat_side <- rbind(latData, radData)
corrDat_PMI <- dcast(dat_side, PPT+GRP+DOM ~ TASK, value.var = 'PMI')
dat <- aggregate(PMI ~ TASK * GRP * PPT, mean, data = dat_side) 
corrDat <- dcast(dat, PPT+GRP ~ TASK, value.var = 'PMI')

ggscatter(corrDat_PMI, x = 'Lateral', y = 'Radial', 
          add = 'reg.line', conf.int = FALSE, add.params = list(color = "black"),
          cor.coef = TRUE, size = 1.5, cor.coef.size = 3, cor.method = 'spearman') + 
  facet_wrap(~GRP*DOM) + 
  ylab('Radial reaching PMI (deg)') + xlab('Lateral reaching PMI (deg)')

ggsave('reaching_correlations_spearman.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

##### dPMI data ######
datDir <- rbind(latData_dir, radData_dir)
corrDat_dir <- dcast(dat_side, PPT+GRP+SIDE ~ TASK, value.var = 'PMI')

ggscatter(corrDat_dir, x = 'Lateral', y = 'Radial', 
          add = 'reg.line', conf.int = FALSE, add.params = list(color = "black"),
          cor.coef = TRUE, size = 1.5, cor.coef.size = 3, cor.method = 'spearman') + 
  facet_wrap(~GRP*SIDE)+ 
  ylab('Radial reaching PMI (deg)') + xlab('Lateral reaching PMI (deg)')

ggsave('dir_reaching_correlations_spearman.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)


#file with data organised by participant
colnames(latData)[which(names(latData) == "Free")] <- "LAT_FREE"
colnames(latData)[which(names(latData) == "Peripheral")] <- "LAT_PER"
colnames(latData)[which(names(latData) == "PMI")] <- "LAT_PMI"
#names(latData)[8] <- 'dPMI_lat'
colnames(radData)[which(names(radData) == "Free")] <- "RAD_FREE"
colnames(radData)[which(names(radData) == "Peripheral")] <- "RAD_PER"
colnames(radData)[which(names(radData) == "PMI")] <- "RAD_PMI"
#names(radData)[8] <- 'dPMI_rad'
allDat <- merge(latData, radData, by = c('PPT', 'GRP', 'DOM'), all = TRUE)
colnames(allDat)[which(names(allDat) == "DIAGNOSIS.x")] <- "DIAGNOSIS"

###### TVA data ######
# ranaming to match reaching data-frame
names(TVAData)[1] <- 'PPT'
names(TVAData)[6] <- 'DIAGNOSIS'
# merge with allDat
allDat <- merge(allDat, TVAData, by = c('PPT', 'GRP', 'DIAGNOSIS'), all = TRUE)
corr_allDat <- na.omit(allDat)

setwd(anaPath)
write.csv(corr_allDat, 'all-task_correlate.csv', row.names = FALSE)


##### TVA lateral reaching
### correlate K and lateral reaching
ggscatter(corr_allDat, x = 'LAT_PMI', y = 'K', add = 'reg.line', conf.int = TRUE,
          cor.coef = TRUE, size = 1, cor.coef.size = 3, cor.method = 'spearman') +
  facet_grid(cols = vars(DOM), rows = vars(GRP)) + 
  ylab('vSTM') + xlab('Lateral reaching PMI (deg)') +
  theme(text = element_text(size = 10))

ggsave('LAT-k_side.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

# collapsed across side
KLAT <- aggregate(LAT_PMI ~ K * PPT * GRP, mean, data = corr_allDat)
ggscatter(KLAT, x = 'LAT_PMI', y = 'K', color = 'grey50', 
          add = 'reg.line', conf.int = FALSE, add.params = list(color = "black"),
          cor.coef = TRUE, size = 1.5, cor.coef.size = 3, cor.method = 'spearman') + 
  facet_wrap(~GRP) + 
  ylab('vSTM (# items)') + xlab('Lateral reaching error (deg)') +
  theme(text = element_text(size = 6)) + theme_classic() 

ggsave('LAT-k_grp.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

### correlate C and lateral reaching
ggscatter(corr_allDat, x = 'LAT_PMI', y = 'C', add = 'reg.line', conf.int = TRUE,
          cor.coef = TRUE, size = 1, cor.coef.size = 3, cor.method = 'spearman') + 
  facet_grid(cols = vars(DOM), rows = vars(GRP)) + 
  ylab('Processing speed (items/s)') + xlab('Lateral reaching PMI (deg)') +
  theme(text = element_text(size = 10))

ggsave('LAT-C_side.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

#### EDITED FOR POSTER PRESENTATION
# collapsed across side
CLAT <- aggregate(LAT_PMI ~ C * PPT * GRP, mean, data = corr_allDat)
ggscatter(CLAT, x = 'LAT_PMI', y = 'C', color = 'grey50', 
          add = 'reg.line', conf.int = FALSE, add.params = list(color = "black"),
          cor.coef = TRUE, size = 1.5, cor.coef.size = 3, cor.method = 'spearman') + 
  facet_wrap(~GRP) + 
  ylab('Processing speed (items/s)') + xlab('Lateral reaching error (deg)') +
  theme(text = element_text(size = 6)) + theme_classic() 

ggsave('LAT-C_grp.png', plot = last_plot(), device = NULL, dpi = 300, 
       width = 3.5, height = 3, path = anaPath)

### correlate K and radial reaching
ggscatter(corr_allDat, x = 'PMI_rad', y = 'K', add = 'reg.line', conf.int = TRUE,
          cor.coef = TRUE, size = 1, cor.coef.size = 3, cor.method = 'spearman') + 
  facet_grid(cols = vars(SIDE), rows = vars(GRP)) + 
  ylab('vSTM') + xlab('Radial reaching PMI (deg)') +
  theme(text = element_text(size = 10))

ggsave('RAD-k_side.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

# collapsed across side
KRAD <- aggregate(PMI_rad ~ K * PPT * GRP, mean, data = corr_allDat)
ggscatter(KRAD, x = 'PMI_rad', y = 'K', add = 'reg.line', conf.int = TRUE,
          cor.coef = TRUE, size = 1, cor.coef.size = 3, cor.method = 'spearman') + 
  facet_wrap(~GRP) + 
  ylab('vSTM') + xlab('Radial reaching PMI (deg)') +
  theme(text = element_text(size = 10))

ggsave('RAD-k_grp.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

### correlate C and radial reaching
ggscatter(corr_allDat, x = 'PMI_rad', y = 'C', add = 'reg.line', conf.int = TRUE,
          cor.coef = TRUE, size = 1, cor.coef.size = 3, cor.method = 'spearman') + 
  facet_grid(cols = vars(SIDE), rows = vars(GRP)) + 
  ylab('Processing speed (items/s)') + xlab('Radial reaching PMI (deg)') +
  theme(text = element_text(size = 10))

ggsave('RAD-C_side.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

# collapsed across side
CRAD <- aggregate(RAD_PMI ~ C * PPT * GRP, mean, data = corr_allDat)
ggscatter(CRAD, x = 'RAD_PMI', y = 'C', add = 'reg.line', conf.int = TRUE,
          cor.coef = TRUE, size = 1, cor.coef.size = 3, cor.method = 'spearman') +  
  facet_wrap(~GRP) + 
  ylab('Processing speed (items/s)') + xlab('Radial reaching PMI (deg)') +
  theme(text = element_text(size = 10))

ggsave('RAD-C_grp.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)


