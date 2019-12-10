##########################################################################################################
## ANALYSIS OF TVA FITS DATA MODELLED IN MATLAB ##
## AG.Mitchell 06.12.2019 ##

##### file info #####
#on mac
#dataPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/data/formatted_TVA'
#fitPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/analysis/TVA/fits/'
#anaPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/analysis/TVA/'
# on pc
dataPath <- 'S:/groups/DMT/data/formatted_TVA/'
fitPath <- ("S:/groups/DMT/analysis/TVA/fits/") # Enter path to data
anaPath <- "S:/groups/DMT/analysis/TVA/"

# Enter directory to save converted files to
setwd(fitPath)

##### functions #####
library(stringr)
library(ggplot2)

# list all files in working directory
txt_filelist <- list.files(fitPath, ".txt")

# converting to .csv files
for (i in 1:length(txt_filelist)) {
  FILE = read.table(file=txt_filelist[i], header = TRUE, sep = '\t')
  
  write.table(
    FILE,
    file=paste0(fitPath,
                sub(".txt","",txt_filelist[i]),".csv"),
    row.names=FALSE,
    quote = FALSE,
    sep=",")
}
 
##### extracting key data #####
# create data-frame
tva_values <- read.csv(text = 'SUB,COND,K,C,T0')
tva_filelist <- dir(fitPath, recursive = TRUE, full.names = FALSE, pattern = '.csv')
##### BOTH ECCs #####
# key values from .csv files - ALL COND FIRST
for (file in tva_filelist){
  if (isTRUE(str_sub(basename(file), -10, -10)=='r')){
    tmp <- read.csv(file)[, c(2,15,28)]
    tmp$COND <- 'all'
    tmp$SUB <- substr(basename(file), 8, 10)
    tmp <- tmp[, c(5,4,1:3)]
    tva_values <- rbind(tva_values, tmp)
  }
  
}

tva_values$GRP <- factor(substr(tva_values$SUB, 1, 1))
levels(tva_values$GRP) <- c('Control', 'Patient')
tva_values$SUB <- factor(tva_values$SUB)
names(tva_values)[5] <- 't0'

#save tva-values
setwd(anaPath)
write.csv(tva_values, 'tva_data.csv', row.names = FALSE)

##### plotting #####
# processing speed
ggplot(tva_values, aes(x = GRP, y = C)) + 
  geom_jitter(aes(colour = GRP), position = position_jitter(0.1)) + 
  scale_color_manual(values = c('black', 'grey50')) +
  labs(title = 'Processing speed (C)', x = 'Group', y = 'C (item/s)', 
                              element_text(size = 12)) +
  theme_bw() + theme(legend.position = 'none', text = element_text(size = 12))

ggsave('processing-speed.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

# vSTM
ggplot(tva_values, aes(x = GRP, y = K)) + 
  geom_jitter(aes(colour = GRP), position = position_jitter(0.1)) + 
  scale_color_manual(values = c('black', 'grey50')) +
  labs(title = 'vSTM', x = 'Group', y = 'K (number of items)', 
       element_text(size = 12)) +
  theme_bw() + theme(legend.position = 'none', text = element_text(size = 12))

ggsave('vSTM.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

# t0
ggplot(tva_values, aes(x = GRP, y = t0)) + 
  geom_jitter(aes(colour = GRP), position = position_jitter(0.1)) + 
  scale_color_manual(values = c('black', 'grey50')) +
  labs(title = 'Perceptual threshold (t0)', x = 'Group', y = 't0 (ms)', 
       element_text(size = 12)) +
  theme_bw() + theme(legend.position = 'none', text = element_text(size = 12))

ggsave('perceptual-thresh.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

##### ECC-3 #####
setwd(fitPath)
tva_values_ecc3 <- read.csv(text = 'SUB,COND,K,C,T0')
# key values from .csv files - ECC-3
for (file in tva_filelist){
  if (isTRUE(str_sub(basename(file), -10, -10)=='3')){
    tmp <- read.csv(file)[, c(2,15,28)]
    tmp$COND <- '3'
    tmp$SUB <- substr(basename(file), 8, 10)
    tmp <- tmp[, c(5,4,1:3)]
    
    tva_values_ecc3 <- rbind(tva_values_ecc3, tmp)
  }
  
}

# remove spurious result (C = 1000, no fit)
tva_values_ecc3 <- tva_values_ecc3[tva_values_ecc3$C < 1000 ,]

##### ECC-9 #####
tva_values_ecc9 <- read.csv(text = 'SUB,COND,K,C,T0')
# key values from .csv files - ECC-3
for (file in tva_filelist){
  if (isTRUE(str_sub(basename(file), -10, -10)=='9')){
    tmp <- read.csv(file)[, c(2,15,28)]
    tmp$COND <- '9'
    tmp$SUB <- substr(basename(file), 8, 10)
    tmp <- tmp[, c(5,4,1:3)]
    
    tva_values_ecc9 <- rbind(tva_values_ecc9, tmp)
  }
  
}

#r-bind both
tva_values_all <- rbind(tva_values_ecc3, tva_values_ecc9)

tva_values_all$GRP <- factor(substr(tva_values_all$SUB, 1, 1))
levels(tva_values_all$GRP) <- c('Control', 'Patient')
tva_values_all$SUB <- factor(tva_values_all$SUB)
names(tva_values_all)[5] <- 't0'

#save tva-values
setwd(anaPath)
write.csv(tva_values_all, 'tva_data_ECC.csv', row.names = FALSE)

##### PLOTTING ######
# processing speed
ggplot(tva_values_all, aes(x = COND, y = C, colour = GRP)) + 
  geom_point(aes(colour = GRP), shape = 1, size = 2) + 
  scale_color_manual(values = c('black', 'grey50')) +
  geom_line(aes(group = SUB), size = 0.5, alpha = 0.5) +
  facet_wrap(~GRP) +
  labs(title = 'Processing speed (C)', x = 'Eccentricity (deg)', y = 'C (item/s)', 
       element_text(size = 12)) +
  theme_bw() + theme(legend.position = 'none', text = element_text(size = 12))

ggsave('ECCprocessing-speed.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

# find t0 for each participant (from other data set, don't need for now)
setwd(dataPath)
dat_filelist <- dir(dataPath, recursive = TRUE, full.names = TRUE, pattern = '.dat')
for (file in dat_filelist){
  if (isTRUE(str_sub(basename(file), -5, -5)=="r")){
    tmp <- read.table(file, header=FALSE, sep="\t")
    timing <- tmp[,2] #column of timing variables
    t <- min(timing, na.rm = TRUE) #getting minimum values
  }
}
