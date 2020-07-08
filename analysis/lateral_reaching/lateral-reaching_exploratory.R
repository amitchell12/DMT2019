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
# on desktop mac
#anaPath <- '/Users/Alex/Documents/DMT/analysis/lateral_reaching'
#dataPath <- '/Users/Alex/Documents/DMT/data'
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

###### DIRECTIONAL ERROR: PMI ######
dir_medians <- aggregate(xerr_mm ~ PPT * VIEW * SIDE * POSITION * SITE * GRP * DIAGNOSIS, 
                         median, data = res)
colnames(dir_medians)[colnames(dir_medians)=='xerr_mm'] <- 'xerr_med' #change name to be more logical
dir_means <- aggregate(xerr_med ~ PPT* VIEW * SIDE * SITE * GRP * DIAGNOSIS, 
                       mean, data = dir_medians)
colnames(dir_means)[colnames(dir_means) == 'xerr_med'] <- 'xerr_mean'

# PMI for directional data (DMI - directional misreaching index)
dPMIdata <- dcast(dir_means, PPT+SIDE+DIAGNOSIS+SITE ~ VIEW) #different data-frame
dPMIdata$PMI <- dPMIdata$Peripheral - dPMIdata$Free
write.csv(dPMIdata, 'lateral-reaching_dirPMI.csv', row.names = FALSE)

# plotting this
ggplot(dir_means, aes(x = VIEW, y = xerr_mean, colour = DIAGNOSIS)) +
  geom_point(shape = 1, size = 2) +
  geom_line(aes(group = PPT), size = 0.5, alpha = .5) +
  facet_grid(cols = vars(SIDE), rows = vars(DIAGNOSIS)) + 
  labs(x = '', y = 'Directional error (mm)', element_text(size = 12)) +
  theme_bw() + theme(legend.position = 'none', text = element_text(size = 10),
                     strip.text.x = element_text(size = 8)) 

# dPMI plot 
ggplot(dPMIdata, aes(x = SIDE, y = PMI, colour = DIAGNOSIS)) + 
  geom_point(shape = 1, size = 1.5, stroke = .8) +
  geom_line(aes(group = PPT), alpha = .5, size = .5) +
  stat_summary(aes(y = PMI, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1, size = 3, group = 1) +
  labs(title = 'Lateral Reaching', x = 'Side', y = 'dPMI (mm)', 
                    element_text(size = 12)) +
  facet_wrap(~DIAGNOSIS) +
  theme_bw() + theme(legend.position = 'none', text = element_text(size = 10),
                     strip.text.x = element_text(size = 10))

ggsave('dPMI_plot.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

## summary dPMI
mean_dPMI <- summarySE(dPMIdata, measurevar = 'PMI', groupvar = c('DIAGNOSIS', 'SIDE'),
                       na.rm = TRUE)
mean_dPMI_all <- summarySE(dPMIdata, measurevar = 'PMI', groupvar = c('DIAGNOSIS', 'SIDE'),
                           na.rm = TRUE)

## ANOVA ##
dPMIanova <- dPMIdata[dPMIdata$PPT != 212 & dPMIdata$PPT != 407 ,] #removing participants where we only have 1 data-point

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

##### DIRECTIONAL ERROR: ECCENTRICITY #####
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
ggplot(av_ecc, aes(x = ECC, y = xerr_med, group = DIAGNOSIS, colour = DIAGNOSIS)) +
  geom_point(size = 3, position = position_dodge(width = .4)) +
  geom_errorbar(aes(ymin=xerr_med-ci, ymax=xerr_med+ci), 
                width=.4, position = position_dodge(width = .4)) + 
  geom_line(aes(group = DIAGNOSIS), size = 0.7, position = position_dodge(width = .4)) +
  geom_hline(yintercept = 0) + ylim(-30,30) +
  labs(x = 'Directional error (mm)', y = 'Eccentricity (°)') +
  facet_grid(~VIEW) + theme_classic() +
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12)
        )

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
ggplot(MT_means, aes(x = DOM, y = MT, colour = DIAGNOSIS, group = DIAGNOSIS)) +
  geom_point(shape = 1, size = 2) +
  geom_line(aes(group = PPT), size = 0.5, alpha = .5) +
  facet_grid(cols = vars(VIEW), rows = vars(DIAGNOSIS)) + ylim(0, 1000) +
  labs(title = 'Reach Duration', x = 'Side', 
       y = 'Reach duration (ms)', element_text(size = 12)) +
  theme_bw() + theme(legend.position = 'none', 
                     text = element_text(size = 10),
                     strip.text.x = element_text(size = 10)) 

# average across sides
MTav$DIAGNOSIS <- factor(MTav$DIAGNOSIS, levels = c('HC','MCI','AD'))
ggplot(MTav, aes(x = VIEW, y = MT, colour = DIAGNOSIS, group = PPT)) +
  geom_point(shape = 16, size = 2, position = position_dodge(width = .3)) +
  geom_line(aes(group = PPT), size = 0.5, alpha = .5, 
            position = position_dodge(width = .3)) +
  stat_summary(aes(y = MT, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1, size = 4, group = 1) +
  facet_wrap(~DIAGNOSIS) + 
  labs(title = 'Reach Duration', x = '', 
       y = 'Reach duration (ms)', element_text(size = 12)) +
  theme_classic() + theme(legend.position = 'none', 
                     text = element_text(size = 10),
                     strip.text.x = element_text(size = 10)
                     ) 

ggsave('lateral-MT.png', plot = last_plot(),  device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

# by eccentricity
MT_medians$ECC <- dir_medians$POSITION
#making left side negative
index <- MT_medians$SIDE == 'left'
MT_medians$ECC[index] <- -(MT_medians$ECC[index])
MT_medians$ECC <- factor(MT_medians$ECC)

# summary data
MTecc <- summarySE(MT_medians, measurevar = 'MT', 
                    groupvar = c('DIAGNOSIS','ECC','VIEW','SIDE'), na.rm = TRUE)
MTecc$DIAGNOSIS <- factor(av_ecc$DIAGNOSIS, levels = c('HC','MCI','AD'))

ggplot(MTecc, aes(x = ECC, y = MT, group = DIAGNOSIS, colour = DIAGNOSIS)) +
  geom_point(size = 3, position = position_dodge(width = .4)) +
  geom_errorbar(aes(ymin=MT-ci, ymax=MT+ci), 
                width=.4, position = position_dodge(width = .4)) + 
  geom_line(aes(group = DIAGNOSIS), size = 0.7, position = position_dodge(width = .4)) +
  labs(x = 'Movement time (ms)', y = 'Eccentricity (°)') +
  facet_grid(~VIEW) + theme_classic() +
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12)
  )

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

###### REACTION TIME ######


### analysis of response times
# do in the same manner as with absolute error
# medians
res_offset_medians <- aggregate(
  time_touch_offset ~ ecc * side * task * subject_nr * site * group, median, data = res)

# means of medians
res_offset_means <- aggregate(
  time_touch_offset ~ task * side * subject_nr * site * group, mean, data = res_offset_medians)


# plotting means
# offset
ggplot(res_offset_means, aes(x = side, y = time_touch_offset, colour = site)) +
  geom_point(shape = 1, size = 2) +
  geom_line(aes(group = subject_nr), size = 0.5, alpha = .5) +
  facet_grid(cols = vars(task), rows = vars(diagnosis)) + ylim(0, 1000) +
  labs(title = 'Touch Offset', x = 'Side', 
       y = 'Touch offset RT (ms)', element_text(size = 12)) +
  scale_colour_manual(values = c('grey50', 'black')) +
  theme_bw() + theme(legend.position = 'none', text = element_text(size = 10),
                     strip.text.x = element_text(size = 10)) -> touchoffsetPlot

ggsave('touchoffset_meansPlot.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, width = 5, height = 6.5, path = anaPath)
  
#reach



#means of both sides
res_reach_meansall <- aggregate(reach_duration~task * subject_nr * site * diagnosis, mean, 
                                data= res_reach_means) 
res_reach_meansall$task <- with(res_reach_meansall, factor(task, levels = rev(levels(task))))

ggplot(res_reach_meansall, aes(x = task, y = reach_duration, colour = site)) +
  geom_point(shape = 1, size = 2) +
  geom_line(aes(group = subject_nr), size = 0.5, alpha = .5) +
  facet_wrap(~diagnosis) + ylim(0, 1000) +
  labs(title = 'Lateral reaching', x = 'Task', 
       y = 'Reach duration (ms)', element_text(size = 12)) +
  scale_colour_manual(values = c('grey50', 'black')) +
  theme_bw() + theme(legend.position = 'none', text = element_text(size = 10),
                     strip.text.x = element_text(size = 10)) -> reachPlot

ggsave('reachDur_means.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)


# correlating peripheral reach duration with PMI
# cast 
rt_offset <- dcast(res_offset_means, subject_nr+diagnosis+site+side ~ task)
rt_reach <- dcast(res_reach_means, subject_nr+diagnosis+site+side ~ task) #different data-frame
#correlations
corrData <- data.frame(PMIdata$PMI)
colnames(corrData)[colnames(corrData) == 'PMIdata.PMI'] <- 'PMI' #renaming
corrData$pAE <- PMIdata$periph
corrData$reachRT <- rt_reach$periph
corrData$offsetRT <- rt_offset$periph
corrData$diagnosis <- PMIdata$diagnosis
corrData$side <- PMIdata$side

# reach dur plot
ggscatter(corrData, x = "pAE", y = "reachRT", add = 'reg.line', conf.int = TRUE,
          cor.coef = TRUE, cor.method = 'pearson') + 
  facet_grid(cols = vars(side), rows = vars(diagnosis))
# touch offset plot
ggscatter(corrData, x = "pAE", y = "offsetRT", add = 'reg.line', conf.int = TRUE,
          cor.coef = TRUE, cor.method = 'pearson') + 
  facet_grid(cols = vars(side), rows = vars(diagnosis))


########### next steps: comparing patients to controls
