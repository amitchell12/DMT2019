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
library(reshape2)
library(ggforce)
library(xpose)

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
tva_values <- read.csv(
  text = 'SUB,ECC,K,C,t0,MU,t1,t2,t3,t4,t5,t6,t7'
  )
tva_filelist <- dir(fitPath, recursive = TRUE, full.names = FALSE, pattern = '.csv')
##### BOTH ECCs #####
# key values from .csv files - ALL COND FIRST
for (file in tva_filelist){
  if (isTRUE(str_sub(basename(file), -10, -10)=='r')){
    tmp <- read.csv(file)[, c(2,15,10,13,28:34)]
    tmp$ECC <- 'all'
    tmp$SUB <- substr(basename(file), 8, 10)
    tmp <- tmp[, c(12,13,1:11)]
    tva_values <- rbind(tva_values, tmp)
  }
  
}

tva_values$GRP <- factor(substr(tva_values$SUB, 1, 1))
levels(tva_values$GRP) <- c('Control', 'Patient')
tva_values$SUB <- factor(tva_values$SUB)
# placing group higher up
tva_values <- tva_values[, c(1,2,14,3:13)]

#adding processing speed
#adding mu to ExpDurc6 and c7
tva_values$ExpDurC6 <- tva_values$ExpDurC6 + tva_values$mu
tva_values$ExpDurC7 <- tva_values$ExpDurC7 + tva_values$mu
#melting so have information by condition
tva_dat <- ans <- melt(
  tva_values, 
  id.vars = c('ECC','SUB','GRP','K','C','t0','mu'), 
  measure.vars = c('ExpDurC1','ExpDurC2','ExpDurC3', 'ExpDurC4', 'ExpDurC5', 'ExpDurC6', 'ExpDurC7'))
#sort new melted data frame
tva_dat <- tva_dat[order(tva_dat$SUB) ,]
#ranaming conditions (1-7) and columns, before adding other columns of value
colnames(tva_dat)[colnames(tva_dat)=="variable"] <- "COND"
tva_dat$COND <- factor(substr(tva_dat$COND, 8, 8))
colnames(tva_dat)[colnames(tva_dat)=="value"] <- "DUR"

## other important variables (predicted k) - compile into variable and add to data-frame
#empty t-tmp data-frame
ttmp <- data.frame(MS=numeric(),
                   pMS=numeric(), 
                   SUB=factor(),
                   COND=factor(),
                 stringsAsFactors=FALSE) 

for (file in tva_filelist){
  if (isTRUE(str_sub(basename(file), -10, -10)=='r')){
    MS <- t(read.csv(file)[, c(81:87)])
    pMS <- t(read.csv(file)[, c(88:94)])
    row.names(MS) <- NULL
    row.names(pMS) <- NULL
    tmp <- as.data.frame(MS)
    tmp$pMS <- pMS
    tmp$SUB <- substr(basename(file), 8, 10)
    tmp$COND <- matrix(1:7, nrow = 7, ncol = 1)
    ttmp <- rbind(ttmp,tmp)
    
  }
  
}

# merging the two data-frames by SUB and COND
tva_dat <- merge(tva_dat, ttmp, by = c('SUB','COND'), all.y = TRUE)
# renaming
colnames(tva_dat)[colnames(tva_dat)=='V1'] <- 'MS'


#save tva-values
setwd(anaPath)
write.csv(tva_dat, 'tva_data.csv', row.names = FALSE)

# plotting predicted duration for each participant
# seperate plots for every 4 participants
# 101-104
ggplot(tva_dat, aes(x = DUR, y = MS)) + 
  geom_point(size = 2) + geom_line(aes(x = DUR, y = pMS), size = 0.5) + 
  facet_wrap_paginate(~SUB, nrow = 2, ncol = 2, scales = 'free_x',
                      strip.position = 'top', page = 1) +
  labs(x = 'Perceived Duration (ms)', y = 'vSTM', 
       element_text(size = 6)) + theme_bw() + 
  theme(legend.position = 'none', text = element_text(size = 10))

ggsave('fits1.pdf', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

# 105-108
ggplot(tva_dat, aes(x = DUR, y = MS)) + 
  geom_point(size = 2) + geom_line(aes(x = DUR, y = pMS), size = 0.5) + 
  facet_wrap_paginate(~SUB, nrow = 2, ncol = 2, scales = 'free_x',
                      strip.position = 'top', page = 2) +
  labs(x = 'Perceived Duration (ms)', y = 'vSTM', 
       element_text(size = 6)) + theme_bw() + 
  theme(legend.position = 'none', text = element_text(size = 10))

ggsave('fits2.pdf', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

# 109-112
ggplot(tva_dat, aes(x = DUR, y = MS)) + 
  geom_point(size = 2) + geom_line(aes(x = DUR, y = pMS), size = 0.5) + 
  facet_wrap_paginate(~SUB, nrow = 2, ncol = 2, scales = 'free_x',
                      strip.position = 'top', page = 3) +
  labs(x = 'Perceived Duration (ms)', y = 'vSTM', 
       element_text(size = 6)) + theme_bw() + 
  theme(legend.position = 'none', text = element_text(size = 10))

ggsave('fits3.pdf', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

# 113-116
ggplot(tva_dat, aes(x = DUR, y = MS)) + 
  geom_point(size = 2) + geom_line(aes(x = DUR, y = pMS), size = 0.5) + 
  facet_wrap_paginate(~SUB, nrow = 2, ncol = 2, scales = 'free_x',
                      strip.position = 'top', page = 3) +
  labs(x = 'Perceived Duration (ms)', y = 'vSTM', 
       element_text(size = 6)) + theme_bw() + 
  theme(legend.position = 'none', text = element_text(size = 10))

ggsave('fits4.pdf', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

# 117-120
ggplot(tva_dat, aes(x = DUR, y = MS)) + 
  geom_point(size = 2) + geom_line(aes(x = DUR, y = pMS), size = 0.5) + 
  facet_wrap_paginate(~SUB, nrow = 2, ncol = 2, scales = 'free_x',
                      strip.position = 'top', page = 3) +
  labs(x = 'Perceived Duration (ms)', y = 'vSTM', 
       element_text(size = 6)) + theme_bw() + 
  theme(legend.position = 'none', text = element_text(size = 10))

ggsave('fits5.pdf', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

# 121-124
ggplot(tva_dat, aes(x = DUR, y = MS)) + 
  geom_point(size = 2) + geom_line(aes(x = DUR, y = pMS), size = 0.5) + 
  facet_wrap_paginate(~SUB, nrow = 2, ncol = 2, scales = 'free_x',
                      strip.position = 'top', page = 3) +
  labs(x = 'Perceived Duration (ms)', y = 'vSTM', 
       element_text(size = 6)) + theme_bw() + 
  theme(legend.position = 'none', text = element_text(size = 10))

ggsave('fits6.pdf', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

# 201-204
ggplot(tva_dat, aes(x = DUR, y = MS)) + 
  geom_point(size = 2) + geom_line(aes(x = DUR, y = pMS), size = 0.5) + 
  facet_wrap_paginate(~SUB, nrow = 2, ncol = 2, scales = 'free_x',
                      strip.position = 'top', page = 3) +
  labs(x = 'Perceived Duration (ms)', y = 'vSTM', 
       element_text(size = 6)) + theme_bw() + 
  theme(legend.position = 'none', text = element_text(size = 10))

ggsave('fits7.pdf', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

##### plotting outcome vars #####
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
tva_values_ecc3 <- read.csv(text = 'SUB,COND,K,C,T0,MU')
# key values from .csv files - ECC-3
for (file in tva_filelist){
  if (isTRUE(str_sub(basename(file), -10, -10)=='3')){
    tmp <- read.csv(file)[, c(2,15,10,13)]
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
