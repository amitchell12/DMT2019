##########################################################################################################
## ANALYSIS OF TVA FITS DATA MODELLED IN MATLAB ##
## AG.Mitchell 06.12.2019 ##

##### functions #####
library(stringr)
library(ggplot2)
library(reshape2)
library(ggforce)
library(xpose)
library(plyr)
library(ez)
library(psychReport)
library(ggpubr)

##### file info #####
# on pc
#fitPath <- ("S:/groups/DMT/analysis/TVA/all/") # Enter path to data
#anaPath <- "S:/groups/DMT/analysis/TVA/"
# on mac
#fitPath <- "/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/analysis/TVA/fits/" # Enter path to data
#anaPath <- "/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/analysis/TVA/"
#dataPath <- "/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/data/"
#latPath <- "/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/analysis/lateral_reaching/"
# on desktop mac
fitPath <- "/Users/Alex/Documents/DMT/analysis/TVA/all/" # Enter path to data
anaPath <- "/Users/Alex/Documents/DMT/analysis/TVA/all/"
dataPath <- "/Users/Alex/Documents/DMT/data/"
latPath <- "/Users/Alex/Documents/DMT/analysis/lateral_reaching/" 

# Enter directory to save converted files to
setwd(fitPath)

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
output1 <- read.csv('Output.csv') #data where t0 not fixed
output2 <- read.csv('OutputFixedt0.csv') #data where t0 was fixed (<0 initially)
# bind
output <- rbind(output1,output2)

# extract important values from data-set
tva_values <- output[, c(1,2,10,11,15,16,13,28:34,81:94,19)]
# labelling
tva_values$SUB <- factor(substr(tva_values$ID, 8, 10))
tva_values$GRP <- factor(substr(tva_values$SUB, 1, 1))
tva_values$SITE <- factor(substr(tva_values$SUB, 1, 1))
#Norwich control data = 3, convert to 1 for same group
for (i in 1:length(tva_values$GRP)){
  if (isTRUE(tva_values$GRP[i] == '3')){
    tva_values$GRP[i] = 1
  }
}
for (i in 1:length(tva_values$GRP)){
  if (isTRUE(tva_values$GRP[i] == '4')){
    tva_values$GRP[i] = 2
  }
}
# same for site, UOE= 1, UEA = 2
for (i in 1:length(tva_values$SITE)){
  if (isTRUE(tva_values$SITE[i] == '2')){
    tva_values$SITE[i] = 1
  }
  if (isTRUE(tva_values$SITE[i] == '3')){
    tva_values$SITE[i] = 2
  }
  if (isTRUE(tva_values$SITE[i] == '4')){
    tva_values$SITE[i] = 2
  }
}

#factors!
tva_values$GRP <- factor(tva_values$GRP)
tva_values$SITE <- factor(tva_values$SITE)
levels(tva_values$GRP) <- c('Control', 'Patient')
levels(tva_values$SITE) <- c('UOE', 'UEA')

# adding demographic info to this data-frame
setwd(dataPath) 
control_demos <- read.csv('control_demographics.csv')
patient_demos <- read.csv('patient_demographics.csv')
#extracting ACE data into seperate data-frame
ACEscores <- patient_demos[ ,c(1, 8:13)]
#isolating patient demographic information to bind with control
patient_demos <- patient_demos[, c(1:6)]
demo <- rbind(control_demos, patient_demos)
names(demo)[1] <- 'SUB'

# merge tva_values with demographics
tva_values <- merge(demo, tva_values, by = 'SUB')

# count N in groups and missing data
Ngroup <- count(tva_values, vars = 'diagnosis')
SPD <- count(tva_values, vars = 'SPD')

# filter
# remove data with SPD > 0 (poorly fitted)
tva_values <- tva_values[tva_values$SPD == 0 ,]
# re-count N per group without SPD > 0 (poorly fitted data)
Ngroup_filt <- count(tva_values, vars = 'diagnosis')

## reorganise data-frame
tva_values <- tva_values[, c(1:6,36:37,8:34)]
tva_values <- tva_values[order(tva_values$SUB) ,]

# adding mu to ExpDurc6 and c7 (unmasked conditions)
tva_values$ExpDurC6 <- tva_values$ExpDurC6 + tva_values$mu
tva_values$ExpDurC7 <- tva_values$ExpDurC7 + tva_values$mu
# saving tva-values
setwd(anaPath)
write.csv(tva_values, 'tva_values.csv', row.names = FALSE)

##### DATA-FRAME FOR FITTING ######
# creating data-frame so can have predicted and actual values for each condition
# melting so have information by condition
# first melt by letter duration
tva_dat <- melt(
  tva_values, 
  id.vars = c('SUB','diagnosis','SITE', 'K','C','t0','mu'), 
  measure.vars = c('ExpDurC1','ExpDurC2','ExpDurC3', 'ExpDurC4', 'ExpDurC5', 'ExpDurC6', 'ExpDurC7'))
#ranaming conditions (1-7) and columns, before adding other columns of value
colnames(tva_dat)[colnames(tva_dat)=="variable"] <- "COND"
tva_dat$COND <- factor(substr(tva_dat$COND, 8, 8))
colnames(tva_dat)[colnames(tva_dat)=="value"] <- "DUR"

# then melt by mean score
tva_MS <- melt(
  tva_values, 
  id.vars = c('SUB','diagnosis'), 
  measure.vars = c('MeanScoreC1','MeanScoreC2','MeanScoreC3', 'MeanScoreC4', 
                   'MeanScoreC5', 'MeanScoreC6', 'MeanScoreC7'))
#ranaming conditions (1-7) and columns, before adding other columns of value
colnames(tva_MS)[colnames(tva_MS)=="variable"] <- "COND"
tva_MS$COND <- factor(substr(tva_MS$COND, 11, 11))
colnames(tva_MS)[colnames(tva_MS)=="value"] <- "MS"

# melt by pred mean score
tva_PMS <- melt(
  tva_values, 
  id.vars = c('SUB','diagnosis'), 
  measure.vars = c('PredMeanScoreC1','PredMeanScoreC2','PredMeanScoreC3', 'PredMeanScoreC4', 
                   'PredMeanScoreC5', 'PredMeanScoreC6', 'PredMeanScoreC7'))
#ranaming conditions (1-7) and columns, before adding other columns of value
colnames(tva_PMS)[colnames(tva_PMS)=="variable"] <- "COND"
tva_PMS$COND <- factor(substr(tva_PMS$COND, 15, 15))
colnames(tva_PMS)[colnames(tva_PMS)=="value"] <- "PredMS"

## finally: merge them all into one data frame with 7 conds per PP
tva_dat <- merge(tva_dat, tva_MS, by = c('SUB', 'COND', 'diagnosis'))
tva_dat <- merge(tva_dat, tva_PMS, by = c('SUB', 'COND', 'diagnosis'))

#save tva-values
setwd(anaPath)
write.csv(tva_dat, 'tva_fits.csv', row.names = FALSE)

######### DATA ANALYSIS #########
# run ANOVAs on processing speed & vSTM
# processing speed
tva_values$SUB <- factor(tva_values$SUB)
C_ANOVA <- ezANOVA(
  data = tva_values
  , dv = .(C)
  , wid = .(SUB)
  , between = .(diagnosis)
  , type = 3,
  return_aov = TRUE,
  detailed = TRUE
)

C_ANOVA$ANOVA
C_ANOVA$`Mauchly's Test for Sphericity`
C_ANOVA$`Sphericity Corrections`
aovC <- aovEffectSize(ezObj = C_ANOVA, effectSize = "pes")
aovDispTable(aovC)

# pairwise t-test to identify group differences
Cttest <- pairwise.t.test(tva_values$C, tva_values$diagnosis, p.adj = 'bonf')
print(Cttest)

# vSTM
K_ANOVA <- ezANOVA(
  data = tva_values
  , dv = .(K)
  , wid = .(SUB)
  , between = .(diagnosis)
  , type = 3,
  return_aov = TRUE,
  detailed = TRUE
)

K_ANOVA$ANOVA
K_ANOVA$`Mauchly's Test for Sphericity`
K_ANOVA$`Sphericity Corrections`
aovK <- aovEffectSize(ezObj = K_ANOVA, effectSize = "pes")
aovDispTable(aovK)

# pairwise t-test to identify group differences
Kttest <- pairwise.t.test(tva_values$K, tva_values$diagnosis, p.adj = 'bonf')
print(Kttest)

### correlate with ACE
names(ACEscores)[1] <- 'SUB'
tvaACE <- merge(tva_values, ACEscores, by = 'SUB')
tvaACE$ACEall <- as.numeric(as.character(tvaACE$ACEall))

## processing speed
ggscatter(tvaACE, x = 'ACEall', y = 'C', add = 'reg.line', conf.int = TRUE,
          cor.coef = TRUE, size = 1, cor.coef.size = 3, cor.method = 'spearman') +
  ylab('Items/s (C)') + xlab('ACE score (%)') +
  theme(text = element_text(size = 10))
#ACE attention score
ggscatter(tvaACE, x = 'ACEattention', y = 'C', add = 'reg.line', conf.int = TRUE,
          cor.coef = TRUE, size = 1, cor.coef.size = 3, cor.method = 'spearman') +
  ylab('Items/s (C)') + xlab('ACE score (%)') +
  theme(text = element_text(size = 10))

## visual working memory
ggscatter(tvaACE, x = 'ACEall', y = 'K', add = 'reg.line', conf.int = TRUE,
          cor.coef = TRUE, size = 1, cor.coef.size = 3, cor.method = 'spearman') +
  ylab('VSTM (k)') + xlab('ACE score (%)') +
  theme(text = element_text(size = 10))
#ACE memory score
ggscatter(tvaACE, x = 'ACEmemory', y = 'K', add = 'reg.line', conf.int = TRUE,
          cor.coef = TRUE, size = 1, cor.coef.size = 3, cor.method = 'spearman') +
  ylab('IVSTM (k)') + xlab('ACE score (%)') +
  theme(text = element_text(size = 10))


######### PLOTTING ##########
# plotting predicted duration for each participant
# plot one for example control, MCI and AD - for presentation

## get individual participant data
tvafits_HC <- tva_dat[tva_dat$SUB == '115',]
tvafits_MCI <- tva_dat[tva_dat$SUB == '210',]
tvafits_AD <- tva_dat[tva_dat$SUB == '408',]

# changing levels for plot
ggplot() + 
  geom_point(data = tvafits_HC, aes(x = DUR, y = MS), size = 4, colour = 'grey50') + 
  geom_line(data = tvafits_HC, aes(x = DUR, y = PredMS), size = 2, colour = 'grey50') + 
  geom_point(data = tvafits_MCI, aes(x = DUR, y = MS), size = 4, colour = 'goldenrod2') + 
  geom_line(data = tvafits_MCI, aes(x = DUR, y = PredMS), size = 2, colour = 'goldenrod2') + 
  geom_point(data = tvafits_AD, aes(x = DUR, y = MS), size = 4, colour = 'dodgerblue3') + 
  geom_line(data = tvafits_AD, aes(x = DUR, y = PredMS), size = 2, colour = 'dodgerblue3') + 
  labs(x = 'Perceived Duration (ms)', y = 'V-STM') +
  theme_classic() +
  theme(legend.position = 'none', axis.text = element_text(size = 20),
        axis.title = element_text(size = 22))

ggsave('example_tvafits.png', plot = last_plot(), device = NULL, dpi = 300, 
       width = 6, height = 5, path = anaPath)


##### plotting outcome vars #####
# loading case-controls to identify individuals with reaching deficit
setwd(latPath)
case_control <- read.csv('lateral-reaching_case-control.csv')
case_control <- aggregate(DEFICIT ~ PPT*DIAGNOSIS, sum, data = case_control)
names(case_control)[1] <- 'SUB'
names(case_control)[2] <- 'diagnosis'

tva_controls <- tva_values[tva_values$diagnosis == 'HC' ,]
tva_controls$DEFICIT <- 0
tva_cases <- merge(tva_values, case_control, by = c('SUB','diagnosis'))
tva_values <- rbind(tva_controls, tva_cases)

tva_values$diagnosis <- factor(tva_values$diagnosis, levels = c('HC', 'MCI', 'AD'))
tva_values$DEFICIT <- factor(tva_values$DEFICIT)

# processing speed
ggplot(tva_values, aes(x = diagnosis, y = C)) + 
  geom_jitter(aes(colour = diagnosis, shape = DEFICIT),
              position = position_jitter(0.2), size = 2) + 
  scale_color_manual(values = c('grey50', 'goldenrod2', 'dodgerblue3')) +
  scale_shape_manual(values = c(16,1,1)) +
  stat_summary(aes(y = C, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1, size = 3, group = 1) +
  labs(title = 'Processing capacity (C)', x = '', y = 'C (item/s)') +
  theme_classic() + 
  theme(legend.position = 'none', 
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        title = element_text(size = 10))

ggsave('processing-speed.png', plot = last_plot(), device = NULL, dpi = 300, 
       width = 4, height = 4.5, path = anaPath)

# vSTM
ggplot(tva_values, aes(x = diagnosis, y = K)) + 
  geom_jitter(aes(colour = diagnosis, shape = DEFICIT), 
              position = position_jitter(0.2),  size = 3) + 
  stat_summary(aes(y = K, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1, size = 3, group = 1) +
  scale_color_manual(values = c('grey50', 'goldenrod2', 'dodgerblue3')) +
  scale_shape_manual(values = c(16,1,1)) +
  labs(title = 'Visual STM (k)', x = '', y = 'K (# items)') +
  theme_classic() + theme(legend.position = 'none', 
                          axis.text = element_text(size = 10),
                          axis.title = element_text(size = 10),
                          title = element_text(size = 10))

ggsave('vSTM.png', plot = last_plot(), device = NULL, dpi = 300, 
       width = 4, height = 4.5, path = anaPath)

# t0
ggplot(tva_values, aes(x = diagnosis, y = t0)) + 
  geom_jitter(aes(colour = diagnosis, shape = DEFICIT), 
              position = position_jitter(0.2),  size = 3) + 
  stat_summary(aes(y = t0, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1, size = 3, group = 1) +
  scale_color_manual(values = c('grey50', 'goldenrod2', 'dodgerblue3')) +
  scale_shape_manual(values = c(16,1,1)) +
  labs(title = 'Perceptual Threshold', x = '', y = 'Speed (ms)') +
  theme_classic() + theme(legend.position = 'none', 
                          axis.text = element_text(size = 10),
                          axis.title = element_text(size = 10),
                          title = element_text(size = 10))
ggsave('T0.png', plot = last_plot(), device = NULL, dpi = 300, 
         width = 4, height = 4.5, path = anaPath)
