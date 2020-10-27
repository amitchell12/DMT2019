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
anaPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/analysis/lateral_reaching'
dataPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/data/'
# on desktop mac
#anaPath <- '/Users/Alex/Documents/DMT/analysis/lateral_reaching'
#dataPath <- '/Users/Alex/Documents/DMT/data'
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
# removing trial with huge RT
res <- res[res$RT < 6000 ,]

# changing levels to be more informative
levels(res$VIEW) <- c('Free', 'Peripheral')
res$SITE <- factor(res$SITE, labels = c('UOE','UEA'))
res$DIAGNOSIS <- factor(res$DIAGNOSIS)
# to match radial reaching - change name of 'POSITION' to 'ECC'
colnames(res)[colnames(res)=='POSITION'] <- 'ECC' 


###### DIRECTIONAL ERROR (xAxis) ######
dir_medians <- aggregate(xerr_mm ~ PPT * VIEW * SIDE * ECC * SITE * DIAGNOSIS * AGE, 
                         median, data = res)
## need to average across side to translate values so -ve = close to fixation
# change sign for left-sided ang-err - use dlpyr funct for this
dir_medians <- dir_medians %>%
  mutate(xerr_mm = ifelse(SIDE == 'left', -xerr_mm, xerr_mm))

# now can average across sides because -ve is towards fix for both R and L trials
dirECC <- aggregate(xerr_mm ~ PPT * VIEW * SITE * ECC * DIAGNOSIS * AGE, 
                     mean, data = dir_medians)
dirECC$ECC <- factor(dirECC$ECC)
# average across target location
dir_means <- aggregate(xerr_mm ~ PPT * VIEW * SITE * DIAGNOSIS * AGE, 
                      mean, data = dirECC)
dir_means$DIAGNOSIS <- factor(dir_means$DIAGNOSIS, levels = c('HC','MCI','AD'))

## DIRECTIONAL ERROR PLOTS ##
# mean all participants
ggplot(dir_means, aes(x = VIEW, y = xerr_mm, colour = DIAGNOSIS, group = PPT)) +
  geom_hline(yintercept = 0) +
  geom_point(shape = 16, size = 2, position = position_dodge(width = .3)) +
  stat_summary(aes(y = xerr_mm, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1, size = 4, group = 1) +
  geom_line(aes(group = PPT), alpha = .4, size = 0.7, position = position_dodge(width = .3)) +
  facet_grid(~DIAGNOSIS) +
  scale_color_manual(values = c('grey50','grey50','grey50')) +
  labs(x = 'Side', y = 'Directional error (mm)') + 
  theme_classic() + theme(legend.position = 'none',
                          axis.text = element_text(size = 10),
                          axis.title = element_text(size = 12),
                          strip.text = element_text(size = 12)
                          )
  
ggsave('LAT_DIRmeans.png', plot = last_plot(), device = NULL, dpi = 300, 
       width = 5, height = 4, path = anaPath)

## plotting eccentricity
avDIR <- summarySE(dirECC, measurevar = 'xerr_mm', 
                              groupvar = c('DIAGNOSIS','ECC','VIEW'), na.rm = TRUE)
avDIR$DIAGNOSIS <- factor(avDIR$DIAGNOSIS, levels = c('HC','MCI','AD'))
# plot 
ggplot(avDIR, aes(x = ECC, y = xerr_mm, group = DIAGNOSIS, colour = DIAGNOSIS, 
                   shape = DIAGNOSIS)) +
  geom_point(size = 3, position = position_dodge(width = .4)) +
  geom_errorbar(aes(ymin=xerr_mm-ci, ymax=xerr_mm+ci), 
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

ggsave('LAT_DIR_ECC.png', plot = last_plot(), device = NULL, dpi = 300, 
       width = 6, height = 5, path = anaPath)

## ANOVA ##
# FULL ANOVA ON ECCENTRICITY DATA
dirECC$ECC <- factor(dirECC$ECC)
DECC_ANOVA <- ezANOVA(
  data = dirECC
  , dv = .(xerr_mm)
  , wid = .(PPT)
  , within = .(VIEW, ECC)
  , between = .(DIAGNOSIS)
  , between_covariates = .(AGE)
  , type = 3,
  return_aov = TRUE,
  detailed = TRUE
)

DECC_ANOVA$ANOVA
DECC_ANOVA$`Mauchly's Test for Sphericity`
DECC_ANOVA$`Sphericity Corrections`
aovDECC <- aovEffectSize(ezObj = DECC_ANOVA, effectSize = "pes")
aovDispTable(aovDECC)

###### AMPLITUDE ERROR (Y-axis) ######
# as OpenSesame [0,0] is in the centre, change sign for lower targets (height = bottom)
# so when average, all going in same direction
res_y <- res %>%
  mutate(yerr_mm = ifelse(height == 'bottom' | height == 'mid', 
                          -yerr_mm, yerr_mm))

AMP_medians <- aggregate(yerr_mm ~ PPT * VIEW * SIDE * ECC * SITE * DIAGNOSIS * AGE, 
                         median, data = res_y)


# now can average across sides because -ve is towards fix for both R and L trials
AMP_ECC <- aggregate(yerr_mm ~ PPT * VIEW * SITE * ECC * DIAGNOSIS * AGE, 
                    mean, data = AMP_medians)
AMP_ECC$ECC <- factor(AMP_ECC$ECC)
# average across target location
AMP_means <- aggregate(yerr_mm ~ PPT * VIEW * SITE * DIAGNOSIS * AGE, 
                       mean, data = AMP_ECC)
AMP_means$DIAGNOSIS <- factor(AMP_means$DIAGNOSIS, levels = c('HC','MCI','AD'))

## DIRECTIONAL ERROR PLOTS ##
# mean all participants
ggplot(AMP_means, aes(x = VIEW, y = yerr_mm, colour = DIAGNOSIS, group = PPT)) +
  geom_hline(yintercept = 0) +
  geom_point(shape = 16, size = 2, position = position_dodge(width = .3)) +
  stat_summary(aes(y = yerr_mm, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1, size = 4, group = 1) +
  geom_line(aes(group = PPT), alpha = .4, size = 0.7, position = position_dodge(width = .3)) +
  facet_grid(~DIAGNOSIS) +
  scale_color_manual(values = c('grey50','grey50','grey50')) +
  labs(x = 'Side', y = 'Amplitude error (mm)') + 
  theme_classic() + theme(legend.position = 'none',
                          axis.text = element_text(size = 10),
                          axis.title = element_text(size = 12),
                          strip.text = element_text(size = 12)
  )

ggsave('LAT_AMPmeans.png', plot = last_plot(), device = NULL, dpi = 300, 
       width = 5, height = 6, path = anaPath)

## plotting eccentricity
avAMP <- summarySE(AMP_medians, measurevar = 'yerr_mm', 
                    groupvar = c('DIAGNOSIS','ECC','VIEW'), na.rm = TRUE)
avAMP$DIAGNOSIS <- factor(avAMP$DIAGNOSIS, levels = c('HC','MCI','AD'))
avAMP$ECC <- factor(avAMP$ECC)
# plot 
ggplot(avAMP, aes(x = ECC, y = yerr_mm, group = DIAGNOSIS, colour = DIAGNOSIS, 
                   shape = DIAGNOSIS)) +
  geom_point(size = 3, position = position_dodge(width = .4)) +
  geom_errorbar(aes(ymin=yerr_mm-ci, ymax=yerr_mm+ci), 
                width=.4, position = position_dodge(width = .4)) + 
  geom_line(aes(group = DIAGNOSIS), size = 0.7, position = position_dodge(width = .4)) +
  geom_hline(yintercept = 0) + 
  scale_color_manual(values = c('black','grey30','grey60')) +
  labs(x = 'Eccentricity (°)', y = 'Lateral reaching error (y-axis, mm)') +
  facet_grid(~VIEW) + theme_classic() +
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12)
  )

ggsave('LAT_AMP_ECC.png', plot = last_plot(), device = NULL, dpi = 300, 
       width = 6, height = 5, path = anaPath)

## ANOVA ##
# FULL ANOVA ON ECCENTRICITY DATA
AMP_ANOVA <- ezANOVA(
  data = AMP_ECC
  , dv = .(yerr_mm)
  , wid = .(PPT)
  , within = .(VIEW, ECC)
  , between = .(DIAGNOSIS)
  , between_covariates = .(AGE)
  , type = 3,
  return_aov = TRUE,
  detailed = TRUE
)

AMP_ANOVA$ANOVA
AMP_ANOVA$`Mauchly's Test for Sphericity`
AMP_ANOVA$`Sphericity Corrections`
aovAMP <- aovEffectSize(ezObj = AMP_ANOVA, effectSize = "pes")
aovDispTable(aovAMP)


###### MOVEMENT TIME #######
MT_medians <- aggregate(MT ~ PPT * VIEW * SIDE * DOM * ECC * SITE * DIAGNOSIS * AGE, 
                        median, data = res)
# average across side
MT_ECC <- aggregate(MT ~ PPT * VIEW * SITE * ECC * DIAGNOSIS * AGE, 
                     mean, data = MT_medians)
# mean
MT_means <- aggregate(MT ~ PPT * VIEW * SITE * DIAGNOSIS,
                      mean, data = MT_ECC)

## plotting :)
# mean all PP
MT_means$DIAGNOSIS <- factor(MT_means$DIAGNOSIS, levels = c('HC','MCI','AD'))

ggplot(MT_means, aes(x = VIEW, y = MT, colour = SITE, group = PPT)) +
  geom_point(shape = 16, size = 2, position = position_dodge(width = .3)) +
  geom_line(aes(group = PPT), size = 0.5, alpha = .5, 
            position = position_dodge(width = .3)) +
  scale_color_manual(values = c('grey50','dodgerblue3')) +
  stat_summary(aes(y = MT, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1, size = 4, group = 1) +
  facet_wrap(~DIAGNOSIS) + 
  labs(x = '', y = 'Movement time (ms)', element_text(size = 12)) +
  theme_classic() + theme(legend.position = 'none', 
                     axis.text = element_text(size = 10),
                     axis.title = element_text(size = 12),
                     strip.text = element_text(size = 12)
                     ) 

ggsave('LAT_MTmeans.png', plot = last_plot(),  device = NULL, dpi = 300, 
       width = 6, height = 4, path = anaPath)

# by eccentricity
# summary data
MTsum <- summarySE(MT_ECC, measurevar = 'MT', 
                    groupvar = c('DIAGNOSIS','ECC','VIEW'), na.rm = TRUE)
MTsum$DIAGNOSIS <- factor(MTsum$DIAGNOSIS, levels = c('HC','MCI','AD'))
MTsum$ECC <- factor(MTsum$ECC)

ggplot(MTsum, aes(x = ECC, y = MT, group = DIAGNOSIS, colour = DIAGNOSIS,
                  shape = DIAGNOSIS)) +
  geom_point(size = 3, position = position_dodge(width = .4)) +
  geom_errorbar(aes(ymin=MT-ci, ymax=MT+ci), 
                width=.4, position = position_dodge(width = .4)) + 
  geom_line(aes(group = DIAGNOSIS), size = 0.7, position = position_dodge(width = .4)) +
  labs(x = 'Eccentricity (°)', y = 'Movement time (ms)') +
  scale_color_manual(values = c('black','grey30','grey60')) +
  ylim(300,800) +
  facet_wrap(~VIEW) + theme_classic() +
  theme(legend.position = 'none',
        legend.title = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 10)
  ) -> MTplot
MTplot

ggsave('LAT_MT_ECC.png', plot = last_plot(),  device = NULL, dpi = 300, 
       width = 7.5, height = 5, path = anaPath)


## ANOVA ## 
MT_ECC$ECC <- factor(MT_ECC$ECC)
MT_ANOVA <- ezANOVA(
  data = MT_ECC
  , dv = .(MT)
  , wid = .(PPT)
  , within = .(VIEW, ECC)
  , between = .(DIAGNOSIS)
  , between_covariates = .(AGE)
  , type = 3,
  return_aov = TRUE,
  detailed = TRUE
)

MT_ANOVA$ANOVA
MT_ANOVA$`Mauchly's Test for Sphericity`
MT_ANOVA$`Sphericity Corrections`
aovMT <- aovEffectSize(ezObj = MT_ANOVA, effectSize = "pes")
aovDispTable(aovMT)

#pair-wise t-test
MTttest <- pairwise.t.test(MT_medians$MT, MT_medians$DIAGNOSIS, p.adj = 'bonf')
print(MTttest)

###### REACTION TIME ######
RT_medians <- aggregate(RT ~ PPT * VIEW * SIDE * DOM * ECC * SITE * DIAGNOSIS * AGE, 
                        median, data = res)
# average across side
RT_ECC <- aggregate(RT ~ PPT * VIEW * ECC * SITE * DIAGNOSIS * AGE, 
                        mean, data = RT_medians)
RT_ECC$ECC <- factor(RT_ECC$ECC)
RT_means <- aggregate(RT ~ PPT * VIEW * SITE * DIAGNOSIS,
                      mean, data = RT_ECC)

## plotting :)
# mean all PP
RT_means$DIAGNOSIS <- factor(RT_means$DIAGNOSIS, levels = c('HC','MCI','AD'))

ggplot(RT_means, aes(x = VIEW, y = RT, colour = SITE, group = PPT)) +
  geom_point(shape = 16, size = 2, position = position_dodge(width = .3)) +
  geom_line(aes(group = PPT), size = 0.5, alpha = .5, 
            position = position_dodge(width = .3)) +
  scale_color_manual(values = c('grey50','dodgerblue3')) +
  stat_summary(aes(y = RT, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1, size = 4, group = 1) +
  facet_wrap(~DIAGNOSIS) + 
  labs(x = '', y = 'Reaction time (ms)', element_text(size = 12)) +
  theme_classic() + theme(legend.position = 'none', 
                          axis.text = element_text(size = 10),
                          axis.title = element_text(size = 12),
                          strip.text = element_text(size = 12)
  ) 

ggsave('LAT_RTmean.png', plot = last_plot(),  device = NULL, dpi = 300, 
       width = 6, height = 4, path = anaPath)

# by eccentricity
# summary data
RTsum <- summarySE(RT_ECC, measurevar = 'RT', 
                   groupvar = c('DIAGNOSIS','ECC','VIEW'), na.rm = TRUE)
RTsum$DIAGNOSIS <- factor(RTsum$DIAGNOSIS, levels = c('HC','MCI','AD'))

ggplot(RTsum, aes(x = ECC, y = RT, group = DIAGNOSIS, colour = DIAGNOSIS,
                  shape = DIAGNOSIS)) +
  geom_point(size = 3, position = position_dodge(width = .4)) +
  geom_errorbar(aes(ymin=RT-ci, ymax=RT+ci), 
                width=.4, position = position_dodge(width = .4)) + 
  geom_line(aes(group = DIAGNOSIS), size = 0.7, position = position_dodge(width = .4)) +
  labs(x = 'Eccentricity (°)', y = 'Reaction time (ms)') +
  scale_color_manual(values = c('black','grey30','grey60')) +
  facet_grid(~VIEW) + theme_classic() +
  ylim(300,800) +
  theme(legend.position = c(.12,.85),
        legend.title = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 10)
  ) -> RTplot
RTplot

ggsave('LAT_RT_ECC.png', plot = last_plot(),  device = NULL, dpi = 300, 
        width = 7.5, height = 5, path = anaPath)

## density plot for RT ##
ggplot(res) +
  geom_density(aes(x = RT, colour = SITE, fill = SITE), alpha = .3) +
  facet_wrap(~DIAGNOSIS) +
  theme_classic() + theme(legend.position = 'bottom')

ggsave('LAT_RT_density.png', plot = last_plot(),  device = NULL, dpi = 300, 
       width = 7.5, height = 5, path = anaPath)

## ANOVA ## 
RT_ANOVA <- ezANOVA(
  data = RT_ECC
  , dv = .(RT)
  , wid = .(PPT)
  , within = .(VIEW, ECC)
  , between = .(DIAGNOSIS)
  , between_covariates = .(AGE)
  , type = 3,
  return_aov = TRUE,
  detailed = TRUE
)

RT_ANOVA$ANOVA
RT_ANOVA$`Mauchly's Test for Sphericity`
RT_ANOVA$`Sphericity Corrections`
aovRT <- aovEffectSize(ezObj = RT_ANOVA, effectSize = "pes")
aovDispTable(aovRT)

### FOR PUBLICATION: combine MT + RT plots ####

TimeFig <- ggarrange(MTplot, RTplot,
                    ncol=2, nrow=1,
                    widths = c(1,1),
                    labels = c('a','b'),
                    hjust = -1)
TimeFig
ggsave('LAT_EXPLOR.png', plot = last_plot(),  device = NULL, dpi = 300, 
       width = 8, height = 4, path = anaPath)

##### CORRELATE PMI + MT, RT ######
# load PMI data
PMIdata <- read.csv('lateralPMI_all.csv')
# average across side
PMIdata <- aggregate(PMI ~ PPT*DIAGNOSIS*SITE*AGE, mean, data = PMIdata)
# cast MT and RT
MT <- dcast(MT_means, PPT+DIAGNOSIS+SITE ~ VIEW)
RT <- dcast(RT_means, PPT+DIAGNOSIS+SITE ~ VIEW)
#renaming for mergings
colnames(MT)[colnames(MT) == 'Free'] <- 'FreeMT' 
colnames(MT)[colnames(MT) == 'Peripheral'] <- 'PeripheralMT' 
colnames(RT)[colnames(RT) == 'Free'] <- 'FreeRT' 
colnames(RT)[colnames(RT) == 'Peripheral'] <- 'PeripheralRT' 

# merging MT with PMI
corrData <- merge(PMIdata, MT, by = c('PPT','DIAGNOSIS','SITE'))
# merging RT with PMI
corrData <- merge(corrData, RT, by = c('PPT','DIAGNOSIS','SITE'))
corrData$DIAGNOSIS <- factor(corrData$DIAGNOSIS, levels = c('HC','MCI','AD'))

## correlate PMI w peripheral MT
ggscatter(corrData, x = "PMI", y = "PeripheralMT", colour = 'grey50',
          add = 'reg.line', conf.int = FALSE, add.params = list(color = "black"),
          cor.coef = TRUE, size = 1.5, cor.coef.size = 3, cor.method = 'spearman')+ 
  facet_wrap(~DIAGNOSIS) +
  labs(x = 'PMI (mm)', y = 'Peripheral movement time (ms)') +
  theme_classic() + 
  theme(text = element_text(size = 10)
        )

ggsave('LAT_MT_PMIcorr.png', plot = last_plot(),  device = NULL, dpi = 300, 
       width = 5, height = 5, path = anaPath)

## correlate PMI w peripheral RT
ggscatter(corrData, x = "PMI", y = "PeripheralRT", colour = 'grey50',
          add = 'reg.line', conf.int = FALSE, add.params = list(color = "black"),
          cor.coef = TRUE, size = 1.5, cor.coef.size = 3, cor.method = 'spearman')+ 
  facet_wrap(~DIAGNOSIS) +
  labs(x = 'PMI (mm)', y = 'Peripheral reaction time (ms)') +
  theme_classic() + 
  theme(text = element_text(size = 10))

ggsave('LAT_RT_PMIcorr.png', plot = last_plot(),  device = NULL, dpi = 300, 
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

## aaaand correlate :)
ggscatter(PMIACE, x = 'ACEall', y = 'PMI', 
          add = 'reg.line', conf.int = FALSE, add.params = list(color = "black"),
          cor.coef = TRUE, size = 1.5, cor.coef.size = 3, cor.method = 'spearman') +
  ylab('PMI (deg)') + xlab('ACE score (%)') +
  theme(text = element_text(size = 10))

ggsave('LAT_PMI_ACEcorr.png', plot = last_plot(),  device = NULL, dpi = 300, 
       width = 5, height = 5, path = anaPath)
