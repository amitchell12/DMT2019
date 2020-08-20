library(readr)
library(ggplot2)
library(Rmisc)
library(gridExtra)
library(plyr)
library(tidyverse)
library(reshape2)
library(Hmisc)
library(ggpubr)
library(ez)
library(psychReport)

###### GETTING DATA #######
#on mac
#anaPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/analysis/lateral_reaching'
# on desktop mac
anaPath <- '/Users/Alex/Documents/DMT/analysis/lateral_reaching'
dataPath <- '/Users/Alex/Documents/DMT/data'
#on pc
#dataPath <- 'S:/groups/DMT/data'
#anaPath <- 'S:/groups/DMT/analysis/lateral_reaching'
setwd(anaPath)

res <- read.csv('lateral-reaching_compiled.csv')
## getting dominant + non-dominant sides, for analysis
names(res)[4] <- 'HAND'
res$HAND <- factor(res$HAND, labels = c('left','right'))
# if hand = side, dominant; else non dominant
res$DOM <- as.numeric(res$SIDE == res$HAND) #1 = dominant, 0 = non-dominant
res$DOM <- factor(res$DOM, labels= c('ND','D'))
# change order so dominance up-front
res <- res[, c(1:8,27,9:26)]

# changing levels to be more informative
res$GRP <- factor(res$GRP, labels = c('Control','Patient'))
levels(res$VIEW) <- c('Free', 'Peripheral')
res$SITE <- factor(res$SITE, labels = c('UOE','UEA'))
res$DIAGNOSIS <- factor(res$DIAGNOSIS)

## find outliers and remove ##
xclude <- read.csv('lateraloutliers.csv')
res <- res[!(res$PPT %in% xclude$PPT), ]

###### DIRECTIONAL ERROR: median, means, PMI ######
dir_medians <- aggregate(xerr_mm ~ PPT * VIEW * SIDE * POSITION * SITE * GRP * DIAGNOSIS, 
                         median, data = res)
colnames(dir_medians)[colnames(dir_medians)=='xerr_mm'] <- 'xerr_med' #change name to be more logical
dir_means <- aggregate(xerr_med ~ PPT* VIEW * SIDE * SITE * GRP * DIAGNOSIS, 
                       mean, data = dir_medians)
colnames(dir_means)[colnames(dir_means) == 'xerr_med'] <- 'xerr_mean'

# PMI for directional data (DMI - directional misreaching index)
dPMIdata <- dcast(dir_means, PPT+SIDE+DIAGNOSIS+GRP+SITE ~ VIEW) #different data-frame
dPMIdata$PMI <- dPMIdata$Peripheral - dPMIdata$Free
write.csv(dPMIdata, 'lateral-reaching_dirPMI.csv', row.names = FALSE)

## DIRECTIONAL ERROR PLOTS ##
# dPMI plot 
dPMIdata$DIAGNOSIS <- factor(dPMIdata$DIAGNOSIS, levels = c('HC','MCI','AD'))
ggplot(dPMIdata, aes(x = SIDE, y = PMI, group = PPT, colour = DIAGNOSIS)) + 
  geom_point(shape = 1, size = 2, stroke = .8, position = position_dodge(width = .2)) +
  geom_line(aes(group = PPT), alpha = .5, size = .5, 
            position = position_dodge(width = .2)) +
  stat_summary(aes(y = PMI, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1, size = 3.5, group = 1) +
  scale_color_manual(values = c('grey50','grey50','grey50')) +
  labs(x = 'Side', y = 'x-axis PMI (mm)', 
       element_text(size = 12)) +
  facet_wrap(~DIAGNOSIS) +
  theme_classic() + theme(legend.position = 'none', 
                          axis.title = element_text(size = 12),
                          axis.text = element_text(size = 10),
                          strip.text = element_text(size = 10) 
                     )

ggsave('LATdPMI_plot.png', plot = last_plot(), device = NULL, dpi = 300, 
       width = 5, height = 6, path = anaPath)

# mean error plot
# summary data first
Err_summary <- summarySE(dir_means, measurevar = 'xerr_mean', 
                         groupvar = c('DIAGNOSIS', 'SIDE', 'VIEW'), na.rm = TRUE)
Err_summary$DIAGNOSIS <- factor(Err_summary$DIAGNOSIS, levels = c('HC','MCI','AD'))

ggplot(Err_summary, aes(x = SIDE, y = xerr_mean, colour = DIAGNOSIS, shape = DIAGNOSIS,
                        group = DIAGNOSIS)) +
  geom_point(size = 4, position = position_dodge(width = .3)) +
  geom_line(aes(group = DIAGNOSIS), size = 0.7, position = position_dodge(width = .3)) +
  geom_errorbar(aes(ymin=xerr_mean-ci, ymax=xerr_mean+ci), 
                width=.3, position = position_dodge(width = .3)) +
  scale_color_manual(values = c('black','grey30','grey60')) +
  labs(x = 'Side', y = 'Lateral reaching error (x-axis, mm)') + 
  facet_wrap(~VIEW) +
  theme_classic() + theme(legend.position = 'bottom', 
                          legend.title = element_blank(),
                          axis.text = element_text(size = 10),
                          axis.title = element_text(size = 12),
                          strip.text = element_text(size = 12))
  
ggsave('LATxerr_means_plot.png', plot = last_plot(), device = NULL, dpi = 300, 
       width = 5, height = 6, path = anaPath)


## summary dPMI
mean_dPMI <- summarySE(dPMIdata, measurevar = 'PMI', groupvar = c('DIAGNOSIS', 'SIDE'),
                       na.rm = TRUE)
mean_dPMI_all <- summarySE(dPMIdata, measurevar = 'PMI', groupvar = c('DIAGNOSIS', 'SIDE'),
                           na.rm = TRUE)

## DIRECTIONAL ERROR ANOVAS ##
dPMIanova <- dPMIdata[dPMIdata$PPT != 212 & dPMIdata$PPT != 407,] #removing participants where we only have 1 data-point

# FULL ANOVA ON FILTERED DATA
DPMI_ANOVA <- ezANOVA(
  data = dPMIanova
  , dv = .(PMI)
  , wid = .(PPT)
  , within = .(SIDE)
  , between = .(DIAGNOSIS)
  , type = 3,
  return_aov = TRUE,
  detailed = TRUE
)

DPMI_ANOVA$ANOVA
DPMI_ANOVA$`Mauchly's Test for Sphericity`
DPMI_ANOVA$`Sphericity Corrections`
aovDPMI <- aovEffectSize(ezObj = DPMI_ANOVA, effectSize = "pes")
aovDispTable(aovDPMI)

## DIRECTIONAL ERROR: ECCENTRICITY ##
# adding eccentricity as a function of side to medians
dir_medians$ECC <- dir_medians$POSITION
#making left side negative
index <- dir_medians$SIDE == 'left'
dir_medians$ECC[index] <- -(dir_medians$ECC[index])
dir_medians$ECC <- factor(dir_medians$ECC)

## plotting data frame
av_ecc <- summarySE(dir_medians, measurevar = 'xerr_med', 
                              groupvar = c('DIAGNOSIS','ECC','VIEW','SIDE'), na.rm = TRUE)
av_ecc$DIAGNOSIS <- factor(av_ecc$DIAGNOSIS, levels = c('HC','MCI','AD'))
# plot 
ggplot(av_ecc, aes(x = ECC, y = xerr_med, group = DIAGNOSIS, colour = DIAGNOSIS, 
                   shape = DIAGNOSIS)) +
  geom_point(size = 3, position = position_dodge(width = .4)) +
  geom_errorbar(aes(ymin=xerr_med-ci, ymax=xerr_med+ci), 
                width=.4, position = position_dodge(width = .4)) + 
  geom_line(aes(group = DIAGNOSIS), size = 0.7, position = position_dodge(width = .4)) +
  geom_hline(yintercept = 0) + ylim(-25,25) +
  scale_color_manual(values = c('black','grey30','grey60')) +
  labs(x = 'Eccentricity (°)', y = 'Lateral reaching error (x-axis, mm)') +
  facet_grid(~VIEW) + theme_classic() +
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12)
        )

ggsave('LATxerr_ECC_plot.png', plot = last_plot(), device = NULL, dpi = 300, 
       width = 6, height = 5, path = anaPath)

## ANOVA ##
dECCanova <- dir_medians[dir_medians$PPT != 212 & dir_medians$PPT != 407 ,]
dECCanova$POSITION <- factor(dECCanova$POSITION)

# FULL ANOVA ON ECCENTRICITY DATA
DECC_ANOVA <- ezANOVA(
  data = dECCanova
  , dv = .(xerr_med)
  , wid = .(PPT)
  , within = .(VIEW, SIDE, POSITION)
  , between = .(DIAGNOSIS)
  , type = 3,
  return_aov = TRUE,
  detailed = TRUE
)

DECC_ANOVA$ANOVA
DECC_ANOVA$`Mauchly's Test for Sphericity`
DECC_ANOVA$`Sphericity Corrections`
aovDECC <- aovEffectSize(ezObj = DECC_ANOVA, effectSize = "pes")
aovDispTable(aovDECC)

###### MOVEMENT TIME #######
MT_medians <- aggregate(MT ~ PPT * VIEW * SIDE * DOM * POSITION * SITE * GRP * DIAGNOSIS, 
                        median, data = res)
MT_means <- aggregate(MT ~ PPT * VIEW * SIDE * DOM * SITE * GRP * DIAGNOSIS,
                      mean, data = MT_medians)
# average across side
MTav <- aggregate(MT ~ PPT * VIEW * SITE * GRP * DIAGNOSIS,
                  mean, data = MT_medians)

## plotting :)
# both sides
MTsummary <- summarySE(MT_means, measurevar = 'MT', 
                       groupvar = c('DIAGNOSIS','VIEW','DOM'), na.rm = TRUE)
MTsummary$DIAGNOSIS <- factor(MTsummary$DIAGNOSIS, levels = c('HC','MCI','AD'))
MTsummary$DOM <- factor(MTsummary$DOM, labels = c('Non-dominant','Dominant'))

ggplot(MTsummary, aes(x = VIEW, y = MT, colour = DIAGNOSIS, group = DIAGNOSIS,
                      shape = DIAGNOSIS)) +
  geom_point(size = 3, position = position_dodge(width = .3)) +
  geom_line(aes(group = DIAGNOSIS), size = 0.5, alpha = .5,
            position = position_dodge(width = .3)) +
  geom_errorbar(aes(ymin=MT-ci, ymax=MT+ci), 
                width=.4, position = position_dodge(width = .3)) +
  scale_color_manual(values = c('black','grey30','grey60')) +
  facet_grid(~DOM) + ylim(300,900) +
  labs(x = '', y = 'Movement time (ms)', element_text(size = 12)) +
  theme_classic() + theme(legend.position = 'bottom', 
                          legend.title = element_blank(),
                          axis.text = element_text(size = 10),
                          axis.title = element_text(size = 12),
                          strip.text = element_text(size = 12)
                     ) 

ggsave('LATMTmean_plot.png', plot = last_plot(),  device = NULL, dpi = 300, 
       width = 4, height = 5, path = anaPath)

# average across sides
MTav$DIAGNOSIS <- factor(MTav$DIAGNOSIS, levels = c('HC','MCI','AD'))
ggplot(MTav, aes(x = VIEW, y = MT, colour = DIAGNOSIS, group = PPT)) +
  geom_point(shape = 16, size = 2, position = position_dodge(width = .3)) +
  geom_line(aes(group = PPT), size = 0.5, alpha = .5, 
            position = position_dodge(width = .3)) +
  scale_color_manual(values = c('grey50','grey50','grey50')) +
  stat_summary(aes(y = MT, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1, size = 4, group = 1) +
  facet_wrap(~DIAGNOSIS) + 
  labs(x = '', y = 'Movement time (ms)', element_text(size = 12)) +
  theme_classic() + theme(legend.position = 'none', 
                     axis.text = element_text(size = 10),
                     axis.title = element_text(size = 12),
                     strip.text = element_text(size = 12)
                     ) 

ggsave('LATMTav_plot.png', plot = last_plot(),  device = NULL, dpi = 300, 
       width = 6, height = 4, path = anaPath)

# by eccentricity
MT_medians$ECC <- dir_medians$POSITION
#making left side negative
index <- MT_medians$SIDE == 'left'
MT_medians$ECC[index] <- -(MT_medians$ECC[index])
MT_medians$ECC <- factor(MT_medians$ECC)

# summary data
MTecc <- summarySE(MT_medians, measurevar = 'MT', 
                    groupvar = c('DIAGNOSIS','ECC','VIEW','SIDE'), na.rm = TRUE)
MTecc$DIAGNOSIS <- factor(MTecc$DIAGNOSIS, levels = c('HC','MCI','AD'))

ggplot(MTecc, aes(x = ECC, y = MT, group = DIAGNOSIS, colour = DIAGNOSIS,
                  shape = DIAGNOSIS)) +
  geom_point(size = 3, position = position_dodge(width = .4)) +
  geom_errorbar(aes(ymin=MT-ci, ymax=MT+ci), 
                width=.4, position = position_dodge(width = .4)) + 
  geom_line(aes(group = DIAGNOSIS), size = 0.7, position = position_dodge(width = .4)) +
  labs(x = 'Eccentricity (°)', y = 'Movement time (ms)') +
  scale_color_manual(values = c('black','grey30','grey60')) +
  facet_wrap(~VIEW) + theme_classic() +
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12)
  )

ggsave('LATMTecc_plot.png', plot = last_plot(),  device = NULL, dpi = 300, 
       width = 7.5, height = 5, path = anaPath)


## ANOVA ## 
MTanova <- MT_means[MT_means$PPT != 212 & MT_means$PPT != 407 ,]
# FULL ANOVA ON MOVEMENT TIME DATA
MT_ANOVA <- ezANOVA(
  data = MTanova
  , dv = .(MT)
  , wid = .(PPT)
  , within = .(VIEW, DOM)
  , between = .(DIAGNOSIS)
  , type = 3,
  return_aov = TRUE,
  detailed = TRUE
)

MT_ANOVA$ANOVA
MT_ANOVA$`Mauchly's Test for Sphericity`
MT_ANOVA$`Sphericity Corrections`
aovMT <- aovEffectSize(ezObj = MT_ANOVA, effectSize = "pes")
aovDispTable(aovMT)

# FULL ANOVA ON MOVEMENT TIME DATA BY ECCENTRICITY
MTECCanova <- MT_medians[MT_medians$PPT != 212 & MT_medians$PPT != 407 ,]
MTECCanova$POSITION <- factor(MTECCanova$POSITION)
MTECC_ANOVA <- ezANOVA(
  data = MTECCanova
  , dv = .(MT)
  , wid = .(PPT)
  , within = .(VIEW, DOM, POSITION)
  , between = .(DIAGNOSIS)
  , type = 3,
  return_aov = TRUE,
  detailed = TRUE
)

MTECC_ANOVA$ANOVA
MTECC_ANOVA$`Mauchly's Test for Sphericity`
MTECC_ANOVA$`Sphericity Corrections`
aovMTECC <- aovEffectSize(ezObj = MTECC_ANOVA, effectSize = "pes")
aovDispTable(aovMTECC)

#pair-wise t-test
MTttest <- pairwise.t.test(MT_medians$MT, MT_medians$DIAGNOSIS, p.adj = 'bonf')
print(MTttest)

###### REACTION TIME ######
RT_medians <- aggregate(RT ~ PPT * VIEW * SIDE * DOM * POSITION * SITE * GRP * DIAGNOSIS, 
                        median, data = res)
RT_means <- aggregate(RT ~ PPT * VIEW * SIDE * DOM * SITE * GRP * DIAGNOSIS,
                      mean, data = RT_medians)
# average across side
RTav <- aggregate(RT ~ PPT * VIEW * SITE * GRP * DIAGNOSIS,
                  mean, data = RT_medians)

## plotting :)
# both sides
RTsummary <- summarySE(RT_means, measurevar = 'RT', 
                       groupvar = c('DIAGNOSIS','VIEW','DOM'), na.rm = TRUE)
RTsummary$DIAGNOSIS <- factor(RTsummary$DIAGNOSIS, levels = c('HC','MCI','AD'))
RTsummary$DOM <- factor(RTsummary$DOM, labels = c('Non-dominant','Dominant'))

ggplot(RTsummary, aes(x = VIEW, y = RT, colour = DIAGNOSIS, group = DIAGNOSIS,
                      shape = DIAGNOSIS)) +
  geom_point(size = 3, position = position_dodge(width = .3)) +
  geom_line(aes(group = DIAGNOSIS), size = 0.5, alpha = .5,
            position = position_dodge(width = .3)) +
  geom_errorbar(aes(ymin=RT-ci, ymax=RT+ci), 
                width=.4, position = position_dodge(width = .3)) +
  scale_color_manual(values = c('black','grey30','grey60')) +
  facet_grid(~DOM) + #ylim(300,900) +
  labs(x = 'Side', y = 'Reaction time (ms)', element_text(size = 12)) +
  theme_classic() + theme(legend.position = 'bottom', 
                          legend.title = element_blank(),
                          axis.text = element_text(size = 10),
                          axis.title = element_text(size = 12),
                          strip.text = element_text(size = 12)
  ) 

ggsave('LATRTmean_plot.png', plot = last_plot(),  device = NULL, dpi = 300, 
       width = 4, height = 5, path = anaPath)

# average across sides
RTav$DIAGNOSIS <- factor(RTav$DIAGNOSIS, levels = c('HC','MCI','AD'))
ggplot(RTav, aes(x = VIEW, y = RT, colour = DIAGNOSIS, group = PPT)) +
  geom_point(shape = 16, size = 2, position = position_dodge(width = .3)) +
  geom_line(aes(group = PPT), size = 0.5, alpha = .5, 
            position = position_dodge(width = .3)) +
  scale_color_manual(values = c('grey50','grey50','grey50')) +
  stat_summary(aes(y = RT, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1, size = 4, group = 1) +
  facet_wrap(~DIAGNOSIS) + 
  labs(x = '', y = 'Movement time (ms)', element_text(size = 12)) +
  theme_classic() + theme(legend.position = 'none', 
                          axis.text = element_text(size = 10),
                          axis.title = element_text(size = 12),
                          strip.text = element_text(size = 12)
  ) 

ggsave('LATRTav_plot.png', plot = last_plot(),  device = NULL, dpi = 300, 
       width = 6, height = 4, path = anaPath)

# by eccentricity
RT_medians$ECC <- dir_medians$POSITION
#making left side negative
index <- RT_medians$SIDE == 'left'
RT_medians$ECC[index] <- -(RT_medians$ECC[index])
RT_medians$ECC <- factor(RT_medians$ECC)

# summary data
RTecc <- summarySE(RT_medians, measurevar = 'RT', 
                   groupvar = c('DIAGNOSIS','ECC','VIEW','SIDE'), na.rm = TRUE)
RTecc$DIAGNOSIS <- factor(RTecc$DIAGNOSIS, levels = c('HC','MCI','AD'))

ggplot(RTecc, aes(x = ECC, y = RT, group = DIAGNOSIS, colour = DIAGNOSIS,
                  shape = DIAGNOSIS)) +
  geom_point(size = 3, position = position_dodge(width = .4)) +
  geom_errorbar(aes(ymin=RT-ci, ymax=RT+ci), 
                width=.4, position = position_dodge(width = .4)) + 
  geom_line(aes(group = DIAGNOSIS), size = 0.7, position = position_dodge(width = .4)) +
  labs(x = 'Eccentricity (°)', y = 'Reaction time (ms)') +
  scale_color_manual(values = c('black','grey30','grey60')) +
  facet_grid(~VIEW) + theme_classic() +
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12)
  )

ggsave('LATRTecc_plot.png', plot = last_plot(),  device = NULL, dpi = 300, 
        width = 7.5, height = 5, path = anaPath)

## ANOVA ## 
RTanova <- RT_means[RT_means$PPT != 212 & RT_means$PPT != 407 ,]
# FULL ANOVA ON MOVEMENT TIME DATA
RT_ANOVA <- ezANOVA(
  data = RTanova
  , dv = .(RT)
  , wid = .(PPT)
  , within = .(VIEW, DOM)
  , between = .(DIAGNOSIS)
  , type = 3,
  return_aov = TRUE,
  detailed = TRUE
)

RT_ANOVA$ANOVA
RT_ANOVA$`Mauchly's Test for Sphericity`
RT_ANOVA$`Sphericity Corrections`
aovRT <- aovEffectSize(ezObj = RT_ANOVA, effectSize = "pes")
aovDispTable(aovRT)

# FULL ANOVA ON MOVEMENT TIME DATA BY ECCENTRICITY
RTECCanova <- RT_medians[RT_medians$PPT != 212 & RT_medians$PPT != 407 ,]
RTECCanova$POSITION <- factor(RTECCanova$POSITION)
RTECC_ANOVA <- ezANOVA(
  data = RTECCanova
  , dv = .(RT)
  , wid = .(PPT)
  , within = .(VIEW, DOM, POSITION)
  , between = .(DIAGNOSIS)
  , type = 3,
  return_aov = TRUE,
  detailed = TRUE
)

RTECC_ANOVA$ANOVA
RTECC_ANOVA$`Mauchly's Test for Sphericity`
RTECC_ANOVA$`Sphericity Corrections`
aovRTECC <- aovEffectSize(ezObj = RTECC_ANOVA, effectSize = "pes")
aovDispTable(aovRTECC)

##### CORRELATE PMI + MT, RT ######
# load PMI data
PMIdata <- read.csv('lateralPMI-filtered.csv')
# cast MT and RT
MT <- dcast(MT_means, PPT+DIAGNOSIS+DOM+SITE ~ VIEW)
RT <- dcast(RT_means, PPT+DIAGNOSIS+DOM+SITE ~ VIEW)
#renaming for mergings
colnames(MT)[colnames(MT) == 'Free'] <- 'FreeMT' 
colnames(MT)[colnames(MT) == 'Peripheral'] <- 'PeripheralMT' 
colnames(RT)[colnames(RT) == 'Free'] <- 'FreeRT' 
colnames(RT)[colnames(RT) == 'Peripheral'] <- 'PeripheralRT' 

# merging MT with PMI
corrData <- merge(PMIdata, MT, by = c('PPT','DIAGNOSIS','SITE','DOM'))
# merging RT with PMI
corrData <- merge(corrData, RT, by = c('PPT','DIAGNOSIS','SITE','DOM'))

## correlate peripheral PMI w peripheral MT
ggscatter(corrData, x = "Peripheral", y = "PeripheralMT", colour = 'grey50',
          add = 'reg.line', conf.int = FALSE, add.params = list(color = "black"),
          cor.coef = TRUE, size = 1.5, cor.coef.size = 3, cor.method = 'spearman')+ 
  facet_grid(cols = vars(GRP), rows = vars(DOM)) +
  labs(x = 'PMI (mm)', y = 'Peripheral movement time (ms)') +
  theme_classic() + 
  theme(text = element_text(size = 10)
        )

ggsave('LATMT-PMI_corr_plot.png', plot = last_plot(),  device = NULL, dpi = 300, 
       width = 5, height = 5, path = anaPath)

## correlate peripheral PMI w peripheral RT
ggscatter(corrData, x = "Peripheral", y = "PeripheralRT", colour = 'grey50',
          add = 'reg.line', conf.int = FALSE, add.params = list(color = "black"),
          cor.coef = TRUE, size = 1.5, cor.coef.size = 3, cor.method = 'spearman')+ 
  facet_grid(cols = vars(GRP), rows = vars(DOM)) +
  labs(x = 'PMI (mm)', y = 'Peripheral reaction time (ms)') +
  theme_classic() + 
  theme(text = element_text(size = 10)
        )

ggsave('LATRT-PMI_corr_plot.png', plot = last_plot(),  device = NULL, dpi = 300, 
       width = 5, height = 5, path = anaPath)

##### CORRELATE ACE #####
setwd(dataPath)
patient_demos <- read.csv('patient_demographics.csv')
#extracting ACE data into seperate data-frame
ACEscores <- patient_demos[ ,c(1, 8:13)]
names(ACEscores)[1] <- 'PPT'
setwd(anaPath)
# merging ACE with PMI
PMIACE <- merge(PMIdata, ACEscores, by = 'PPT')
PMIACE$ACEall <- as.numeric(as.character(PMIACE$ACEall))
PMIACE$ACEvisuospatial <- as.numeric(as.character(PMIACE$ACEvisuospatial))

## aaaand correlate :)
ggscatter(PMIACE, x = 'ACEall', y = 'PMI', add = 'reg.line', conf.int = TRUE,
          cor.coef = TRUE, size = 1, cor.coef.size = 3, cor.method = 'spearman') +
  ylab('PMI (deg)') + xlab('ACE score (%)') +
  facet_wrap(~DOM) +
  theme(text = element_text(size = 10))

# average PMI
PMIACE <- aggregate(PMI ~ PPT+DIAGNOSIS+ACEall, mean, data = PMIACE)
ggscatter(PMIACE, x = 'ACEall', y = 'PMI', 
          add = 'reg.line', conf.int = FALSE, add.params = list(color = "black"),
          cor.coef = TRUE, size = 1.5, cor.coef.size = 3, cor.method = 'spearman') +
  ylab('PMI (deg)') + xlab('ACE score (%)') +
  theme(text = element_text(size = 10))

ggsave('LATPMI-ACEcorr_plot.png', plot = last_plot(),  device = NULL, dpi = 300, 
       width = 5, height = 5, path = anaPath)
