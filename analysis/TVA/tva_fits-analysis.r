##########################################################################################################
## ANALYSIS OF TVA FITS DATA MODELLED IN MATLAB ##
## AG.Mitchell 06.12.2019 ##

##### file info #####
#on mac
#dataPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/data/formatted_TVA'
#fitPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/analysis/TVA/fits/'
#anaPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/analysis/TVA/'
# on pc
fitPath <- ("S:/groups/DMT/analysis/TVA/fits/") # Enter path to data
dataPath <- 'S:/groups/DMT/data/formatted_TVA/'
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
tva_values$SITE <- factor(substr(tva_values$SUB, 1, 1))
#Norwich control data = 3, convert to 1 for same group
for (i in 1:length(tva_values$GRP)){
  if (isTRUE(tva_values$GRP[i] == '3')){
    tva_values$GRP[i] = 1
  }
}
for (i in 1:length(tva_values$SITE)){
  if (isTRUE(tva_values$SITE[i] == '2')){
    tva_values$SITE[i] = 1
  }
}
tva_values$GRP <- factor(tva_values$GRP)
tva_values$SITE <- factor(tva_values$SITE)
levels(tva_values$GRP) <- c('Control', 'Patient')
levels(tva_values$SITE) <- c('UOE', 'UEA')
tva_values$SUB <- factor(tva_values$SUB)
# placing group higher up
tva_values <- tva_values[, c(1,2,14:15,3:13)] ##### check this order with new 'site' column

#adding processing speed
#adding mu to ExpDurc6 and c7
tva_values$ExpDurC6 <- tva_values$ExpDurC6 + tva_values$mu
tva_values$ExpDurC7 <- tva_values$ExpDurC7 + tva_values$mu
# saving tva-values
setwd(anaPath)
write.csv(tva_values, 'tva_values.csv', row.names = FALSE)

#melting so have information by condition
tva_dat <- melt(
  tva_values, 
  id.vars = c('ECC','SUB','GRP','SITE', 'K','C','t0','mu'), 
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

setwd(fitPath)
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
# order by group
tva_dat <- tva_dat[order(tva_dat$GRP) ,]


#save tva-values
setwd(anaPath)
write.csv(tva_dat, 'tva_fits.csv', row.names = FALSE)

# plotting predicted duration for each participant
# seperate plots for every 4 participants
###### redo these so you have document for all participants
# 101-104=6
ggplot(tva_dat, aes(x = DUR, y = MS)) + 
  geom_point(size = 1) + geom_line(aes(x = DUR, y = pMS), size = 0.5) + 
  facet_wrap_paginate(~SUB, nrow = 3, ncol = 2, scales = 'free_x',
                      strip.position = 'top', page = 1) +
  labs(x = 'Perceived Duration (ms)', y = 'vSTM', 
       element_text(size = 6)) + theme_bw() + 
  theme(legend.position = 'none', text = element_text(size = 8))

ggsave('fits1.pdf', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

# 107-112
ggplot(tva_dat, aes(x = DUR, y = MS)) + 
  geom_point(size = 1) + geom_line(aes(x = DUR, y = pMS), size = 0.5) + 
  facet_wrap_paginate(~SUB, nrow = 3, ncol = 2, scales = 'free_x',
                      strip.position = 'top', page = 2) +
  labs(x = 'Perceived Duration (ms)', y = 'vSTM', 
       element_text(size = 6)) + theme_bw() + 
  theme(legend.position = 'none', text = element_text(size = 8))

ggsave('fits2.pdf', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

# 113-118
ggplot(tva_dat, aes(x = DUR, y = MS)) + 
  geom_point(size = 1) + geom_line(aes(x = DUR, y = pMS), size = 0.5) + 
  facet_wrap_paginate(~SUB, nrow = 3, ncol = 2, scales = 'free_x',
                      strip.position = 'top', page = 3) +
  labs(x = 'Perceived Duration (ms)', y = 'vSTM', 
       element_text(size = 6)) + theme_bw() + 
  theme(legend.position = 'none', text = element_text(size = 8))

ggsave('fits3.pdf', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

# 119-124
ggplot(tva_dat, aes(x = DUR, y = MS)) + 
  geom_point(size = 1) + geom_line(aes(x = DUR, y = pMS), size = 0.5) + 
  facet_wrap_paginate(~SUB, nrow = 3, ncol = 2, scales = 'free_x',
                      strip.position = 'top', page = 4) +
  labs(x = 'Perceived Duration (ms)', y = 'vSTM', 
       element_text(size = 6)) + theme_bw() + 
  theme(legend.position = 'none', text = element_text(size = 8))

ggsave('fits4.pdf', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

# 201-207
ggplot(tva_dat, aes(x = DUR, y = MS)) + 
  geom_point(size = 1) + geom_line(aes(x = DUR, y = pMS), size = 0.5) + 
  facet_wrap_paginate(~SUB, nrow = 3, ncol = 2, scales = 'free_x',
                      strip.position = 'top', page = 5) +
  labs(x = 'Perceived Duration (ms)', y = 'vSTM', 
       element_text(size = 6)) + theme_bw() + 
  theme(legend.position = 'none', text = element_text(size = 8))

ggsave('fits5.pdf', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

# 208-302
ggplot(tva_dat, aes(x = DUR, y = MS)) + 
  geom_point(size = 1) + geom_line(aes(x = DUR, y = pMS), size = 0.5) + 
  facet_wrap_paginate(~SUB, nrow = 3, ncol = 2, scales = 'free_x',
                      strip.position = 'top', page = 6) +
  labs(x = 'Perceived Duration (ms)', y = 'vSTM', 
       element_text(size = 6)) + theme_bw() + 
  theme(legend.position = 'none', text = element_text(size = 8))

ggsave('fits6.pdf', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

# 304-309
ggplot(tva_dat, aes(x = DUR, y = MS)) + 
  geom_point(size = 1) + geom_line(aes(x = DUR, y = pMS), size = 0.5) + 
  facet_wrap_paginate(~SUB, nrow = 3, ncol = 2, scales = 'free_x',
                      strip.position = 'top', page = 7) +
  labs(x = 'Perceived Duration (ms)', y = 'vSTM', 
       element_text(size = 6)) + theme_bw() + 
  theme(legend.position = 'none', text = element_text(size = 8))

ggsave('fits7.pdf', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

# 312-316
ggplot(tva_dat, aes(x = DUR, y = MS)) + 
  geom_point(size = 1) + geom_line(aes(x = DUR, y = pMS), size = 0.5) + 
  facet_wrap_paginate(~SUB, nrow = 3, ncol = 2, scales = 'free_x',
                      strip.position = 'top', page = 8) +
  labs(x = 'Perceived Duration (ms)', y = 'vSTM', 
       element_text(size = 6)) + theme_bw() + 
  theme(legend.position = 'none', text = element_text(size = 8))

ggsave('fits8.pdf', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

## average predicted duration for group


##### plotting outcome vars #####
# processing speed
ggplot(tva_values, aes(x = GRP, y = C)) + 
  geom_jitter(aes(colour = SITE), position = position_jitter(0.1)) + 
  scale_color_manual(values = c('grey50', 'red')) +
  stat_summary(aes(y = C, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1, size = 4, group = 1) +
  labs(title = 'Processing speed (C)', x = 'Group', y = 'C (item/s)', 
                              element_text(size = 10)) +
  theme_bw() + theme(legend.position = 'bottom', text = element_text(size = 10))

ggsave('processing-speed.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

# vSTM
ggplot(tva_values, aes(x = GRP, y = K)) + 
  geom_jitter(aes(colour = SITE), position = position_jitter(0.1)) + 
  scale_color_manual(values = c('grey50', 'red')) +
  stat_summary(aes(y = K, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1, size = 4, group = 1) +
  labs(title = 'vSTM', x = 'Group', y = 'K (number of items)', 
       element_text(size = 10)) +
  theme_bw() + theme(legend.position = 'bottom', text = element_text(size = 10))

ggsave('vSTM.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

# t0
ggplot(tva_values, aes(x = GRP, y = t0)) + 
  geom_jitter(aes(colour = SITE), position = position_jitter(0.1)) + 
  scale_color_manual(values = c('grey50', 'red')) +
  stat_summary(aes(y = t0, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1, size = 4, group = 1) +
  labs(title = 'Perceptual threshold (t0)', x = 'Group', y = 't0 (ms)', 
       element_text(size = 10)) +
  theme_bw() + theme(legend.position = 'bottom', text = element_text(size = 10))

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
    tmp <- tmp[, c(6,5,1:4)]
    
    tva_values_ecc3 <- rbind(tva_values_ecc3, tmp)
  }
  
}

# remove spurious result (C = 1000, no fit)
tva_values_ecc3 <- tva_values_ecc3[tva_values_ecc3$C < 1000 ,]

##### ECC-9 #####
tva_values_ecc9 <- read.csv(text = 'SUB,COND,K,C,T0,MU')
# key values from .csv files - ECC-3
for (file in tva_filelist){
  if (isTRUE(str_sub(basename(file), -10, -10)=='9')){
    tmp <- read.csv(file)[, c(2,15,10,13)]
    tmp$COND <- '9'
    tmp$SUB <- substr(basename(file), 8, 10)
    tmp <- tmp[, c(6,5,1:4)]
    
    tva_values_ecc9 <- rbind(tva_values_ecc9, tmp)
  }
  
}

#r-bind both
tva_values_all <- rbind(tva_values_ecc3, tva_values_ecc9)

tva_values_all$GRP <- factor(substr(tva_values_all$SUB, 1, 1))
for (i in 1:length(tva_values_all$GRP)){
  if (isTRUE(tva_values_all$GRP[i] == '3')){
    tva_values_all$GRP[i] = 1
  }
}
tva_values_all$GRP <- factor(tva_values_all$GRP)
tva_values_all
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


