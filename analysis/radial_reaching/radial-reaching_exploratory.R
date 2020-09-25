# EXPLORATORY ANALYSIS FOR RADIAL REACHING TASK 
## AG.Mitchell 22.04.20 ##
library(readr)
library(ggplot2)
library(Rmisc)
library(gridExtra)
library(plyr)
library(dplyr)
library(tidyverse)
library(reshape2)
library(Hmisc)
library(ggpubr)
library(ez)
library(psychReport)

###### GETTING DATA ######
#on mac
#anaPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/analysis/radial_reaching'
#dataPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/data'
#desktop mac
anaPath <- "/Users/Alex/Documents/DMT/analysis/radial_reaching/"
dataPath <- '/Users/Alex/Documents/DMT/data'
#on pc
#dataPath <- 'S:/groups/DMT/data'
#anaPath <- 'S:/groups/DMT/analysis/radial_reaching'
setwd(anaPath)

res <- read.csv('all_radial-reaching_compiled.csv')
res$ECC <- factor(abs(res$POSITION)) #adding eccentricity = absolute target position

##### AE BY ECCENTRICITY ALL TARG LOCS #####
# calculating medians
res_medians <- aggregate(AE ~ PPT*ECC*VIEW*SIDE*DOM*DIAGNOSIS*GRP*SITE*AGE, 
                                                       median, data = res)
# averaging across side
resAE <- aggregate(AE ~ PPT*ECC*VIEW*DIAGNOSIS*SITE*AGE, 
                    mean, data = res_medians)
## CALCULATING ANOVA ##
ALLECC_ANOVA <- ezANOVA(
  data = resAE
  , dv = .(AE)
  , wid = .(PPT)
  , within = .(VIEW, ECC)
  , between = .(DIAGNOSIS)
  , between_covariates = .(SITE, AGE)
  , type = 3,
  return_aov = TRUE,
  detailed = TRUE
)

ALLECC_ANOVA$ANOVA
ALLECC_ANOVA$`Mauchly's Test for Sphericity`
ALLECC_ANOVA$`Sphericity Corrections`
aovECCall <- aovEffectSize(ezObj = ALLECC_ANOVA, effectSize = "pes")
aovDispTable(aovECCall)

## PLOT: eccentricity ##
ECCsummary <- summarySE(resAE, measurevar = 'AE', 
                        groupvar = c('DIAGNOSIS', 'VIEW', 'ECC'), na.rm = TRUE)
ECCsummary$DIAGNOSIS <- factor(ECCsummary$DIAGNOSIS, levels = c('HC','MCI','AD'))

ggplot(ECCsummary, aes(x = ECC, y = AE, group = DIAGNOSIS, colour = DIAGNOSIS, 
                       shape = DIAGNOSIS)) +
  geom_point(size = 3, position = position_dodge(-.2)) +
  geom_errorbar(aes(ymin=AE-ci, ymax=AE+ci), 
                width=.4, position = position_dodge(-.2)) + 
  geom_line(aes(group = DIAGNOSIS), size = 0.7, position = position_dodge(-.2)) +
  scale_color_manual(values = c('black','grey30','grey60')) +
  labs(x = 'Eccentricity (mm)', y = 'Radial reaching error (mm)') +
  facet_grid(~VIEW) + theme_classic() +
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12)
  )

ggsave('RADeccentricity-fig.png', plot = last_plot(), device = NULL, dpi = 300, 
       width = 7, height = 5, path = anaPath)


##### ANGULAR ERROR: median, means, PMI #####
ANGmedians <- aggregate(ANG_ERR ~ PPT * VIEW * SIDE * ECC * SITE * GRP * DIAGNOSIS * AGE * ECC, 
                         median, data = res)
## need to average across side to translate values so -ve = close to fixation
# change sign for left-sided ang-err - use dlpyr funct for this
ANGmedians <- ANGmedians %>%
  mutate(ANG_ERR = ifelse(SIDE == 'Left', -ANG_ERR, ANG_ERR))

# now can average across sides because -ve is towards fix for both R and L trials
ANG_ECC <- aggregate(ANG_ERR ~ PPT * VIEW * SITE * ECC * DIAGNOSIS * AGE, 
                      mean, data = ANGmedians)
# average across target location
ANGmeans <- aggregate(ANG_ERR ~ PPT * VIEW * SITE * DIAGNOSIS * AGE, 
                           mean, data = ANG_ECC)
ANGmeans$DIAGNOSIS <- factor(ANGmeans$DIAGNOSIS, levels = c('HC','MCI','AD'))

### PLOTTING ###
## mean all participants
ggplot(ANGmeans, aes(x = VIEW, y = ANG_ERR, group = PPT, colour = DIAGNOSIS)) +
  geom_hline(yintercept = 0) +
  geom_point(size = 3, shape = 1, position = position_dodge(width = .2)) +
  stat_summary(aes(y = ANG_ERR, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1, size = 4, group = 1) +
  geom_line(aes(group = PPT), alpha = .4, size = 0.7, position = position_dodge(width = .2)) +
  scale_color_manual(values = c('grey50','grey50','grey50')) +
  facet_grid(~DIAGNOSIS) +
  labs(x = 'Side', y = 'Angular error (°)') + 
  theme_classic() + theme(legend.position = 'none',
                          axis.text = element_text(size = 10),
                          axis.title = element_text(size = 12),
                          strip.text = element_text(size = 12)
                          )

ggsave('RAD_ANGERRmeans.png', plot = last_plot(), device = NULL, dpi = 300, 
       width = 5, height = 7, path = anaPath)

## plotting data frame
ANGECC_summary <- summarySE(ANG_ECC, measurevar = 'ANG_ERR', 
                         groupvar = c('DIAGNOSIS', 'VIEW', 'ECC'), na.rm = TRUE)
ANGECC_summary$DIAGNOSIS <- factor(ANGECC_summary$DIAGNOSIS, levels = c('HC','MCI','AD'))

ggplot(ANGECC_summary, aes(x = ECC, y = ANG_ERR, colour = DIAGNOSIS, shape = DIAGNOSIS,
                        group = DIAGNOSIS)) +
  geom_hline(yintercept = 0) +
  geom_point(size = 4, position = position_dodge(width = .3)) +
  geom_line(aes(group = DIAGNOSIS), size = 0.7, position = position_dodge(width = .3)) +
  geom_errorbar(aes(ymin=ANG_ERR-ci, ymax=ANG_ERR+ci), 
                width=.3, position = position_dodge(width = .3)) +
  scale_color_manual(values = c('black','grey30','grey60')) +
  labs(x = 'Side', y = 'Angular error (°)') + 
  facet_wrap(~VIEW) +
  theme_classic() + theme(legend.position = 'bottom', 
                          legend.title = element_blank(),
                          axis.text = element_text(size = 10),
                          axis.title = element_text(size = 12),
                          strip.text = element_text(size = 12))

ggsave('RAD_ANGERR_ECC.png', plot = last_plot(), device = NULL, dpi = 300, 
       width = 5, height = 6, path = anaPath)

### ANOVA ###
# FULL ANOVA ON ECCENTRICITY DATA
ANG_ANOVA <- ezANOVA(
  data = ANG_ECC
  , dv = .(ANG_ERR)
  , wid = .(PPT)
  , within = .(VIEW, ECC)
  , between = .(DIAGNOSIS)
  , between_covariates = .(SITE, AGE)
  , type = 3,
  return_aov = TRUE,
  detailed = TRUE
)

ANG_ANOVA$ANOVA
ANG_ANOVA$`Mauchly's Test for Sphericity`
ANG_ANOVA$`Sphericity Corrections`
aovANGECC <- aovEffectSize(ezObj = ANG_ANOVA, effectSize = "pes")
aovDispTable(aovANGECC)

###### AMPLITUDE ERROR ######
AMPmedians <- aggregate(AMP_ERR ~ PPT * VIEW * SIDE * ECC * SITE * GRP * DIAGNOSIS * AGE * ECC, 
                        median, data = res)
# average across side (-ve = closer to PP, +ve = away from PP)
AMP_ECC <- aggregate(AMP_ERR ~ PPT * VIEW * SITE * ECC * DIAGNOSIS * AGE, 
                     mean, data = AMPmedians)
# average across target location
AMPmeans <- aggregate(AMP_ERR ~ PPT * VIEW * SITE * DIAGNOSIS * AGE, 
                      mean, data = AMP_ECC)
AMPmeans$DIAGNOSIS <- factor(AMPmeans$DIAGNOSIS, levels = c('HC','MCI','AD'))

### PLOTTING ###
## mean all participants
ggplot(AMPmeans, aes(x = VIEW, y = AMP_ERR, group = PPT, colour = DIAGNOSIS)) +
  geom_hline(yintercept = 0) +
  geom_point(size = 3, shape = 1, position = position_dodge(width = .2)) +
  stat_summary(aes(y = AMP_ERR, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1, size = 4, group = 1) +
  geom_line(aes(group = PPT), alpha = .4, size = 0.7, position = position_dodge(width = .2)) +
  scale_color_manual(values = c('grey50','grey50','grey50')) +
  facet_grid(~DIAGNOSIS) +
  labs(x = 'Side', y = 'Amplutude error (mm)') + 
  theme_classic() + theme(legend.position = 'none',
                          axis.text = element_text(size = 10),
                          axis.title = element_text(size = 12),
                          strip.text = element_text(size = 12)
  )

ggsave('RAD_AMPERRmeans.png', plot = last_plot(), device = NULL, dpi = 300, 
       width = 5, height = 7, path = anaPath)

## plotting data frame
AMPECC_summary <- summarySE(AMP_ECC, measurevar = 'AMP_ERR', 
                            groupvar = c('DIAGNOSIS', 'VIEW', 'ECC'), na.rm = TRUE)
AMPECC_summary$DIAGNOSIS <- factor(AMPECC_summary$DIAGNOSIS, levels = c('HC','MCI','AD'))

ggplot(AMPECC_summary, aes(x = ECC, y = AMP_ERR, colour = DIAGNOSIS, shape = DIAGNOSIS,
                           group = DIAGNOSIS)) +
  geom_hline(yintercept = 0) +
  geom_point(size = 4, position = position_dodge(width = .3)) +
  geom_line(aes(group = DIAGNOSIS), size = 0.7, position = position_dodge(width = .3)) +
  geom_errorbar(aes(ymin=AMP_ERR-ci, ymax=AMP_ERR+ci), 
                width=.3, position = position_dodge(width = .3)) +
  scale_color_manual(values = c('black','grey30','grey60')) +
  labs(x = 'Side', y = 'Amplitude error (mm)') + 
  facet_wrap(~VIEW) +
  theme_classic() + theme(legend.position = 'bottom', 
                          legend.title = element_blank(),
                          axis.text = element_text(size = 10),
                          axis.title = element_text(size = 12),
                          strip.text = element_text(size = 12))

ggsave('RAD_AMPERR_ECC.png', plot = last_plot(), device = NULL, dpi = 300, 
       width = 5, height = 6, path = anaPath)

### ANOVA ###
# FULL ANOVA ON ECCENTRICITY DATA
AMP_ANOVA <- ezANOVA(
  data = AMP_ECC
  , dv = .(AMP_ERR)
  , wid = .(PPT)
  , within = .(VIEW, ECC)
  , between = .(DIAGNOSIS)
  , between_covariates = .(SITE, AGE)
  , type = 3,
  return_aov = TRUE,
  detailed = TRUE
)

AMP_ANOVA$ANOVA
AMP_ANOVA$`Mauchly's Test for Sphericity`
AMP_ANOVA$`Sphericity Corrections`
aovAMPECC <- aovEffectSize(ezObj = AMP_ANOVA, effectSize = "pes")
aovDispTable(aovAMPECC)

###### MOVEMENT TIME #######
MT_medians <- aggregate(MT ~ PPT * VIEW * SIDE * DOM * ECC * SITE * DIAGNOSIS, 
                        median, data = res)
MT_medians <- aggregate(MT ~ PPT * VIEW * DOM * ECC * SITE * DIAGNOSIS, 
                        median, data = MT_medians)
MT_means <- aggregate(MT ~ PPT * VIEW * DOM * SITE * DIAGNOSIS,
                      mean, data = MT_ECC)
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
  facet_grid(~DOM) + ylim(500,1000) +
  labs(x = '', y = 'Movement time (ms)', element_text(size = 12)) +
  theme_classic() + theme(legend.position = 'bottom', 
                          legend.title = element_blank(),
                          axis.text = element_text(size = 10),
                          axis.title = element_text(size = 12),
                          strip.text = element_text(size = 12)
  ) 

ggsave('RADMTmean_plot.png', plot = last_plot(),  device = NULL, dpi = 300, 
       width = 4, height = 5, path = anaPath)

# average across sides
MTav$DIAGNOSIS <- factor(MTav$DIAGNOSIS, levels = c('HC','MCI','AD'))
ggplot(MTav, aes(x = VIEW, y = MT, colour = SITE, group = PPT)) +
  geom_point(shape = 16, size = 2, position = position_dodge(width = .3)) +
  geom_line(aes(group = PPT), size = 0.5, alpha = .5, 
            position = position_dodge(width = .3)) +
  scale_color_manual(values = c('dodgerblue3','grey50')) +
  stat_summary(aes(y = MT, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1, size = 4, group = 1) +
  facet_wrap(~DIAGNOSIS) + 
  labs(x = '', y = 'Movement time (ms)', element_text(size = 12)) +
  theme_classic() + theme(legend.position = 'bottom', 
                          axis.text = element_text(size = 10),
                          axis.title = element_text(size = 12),
                          strip.text = element_text(size = 12)
  ) 

ggsave('RADMTav_plot.png', plot = last_plot(),  device = NULL, dpi = 300, 
       width = 6, height = 4, path = anaPath)

# by eccentricity
MT_medians$POSITION <- factor(MT_medians$POSITION)

# summary data
MTecc <- summarySE(MT_medians, measurevar = 'MT', 
                   groupvar = c('DIAGNOSIS','POSITION','VIEW'), na.rm = TRUE)
MTecc$DIAGNOSIS <- factor(MTecc$DIAGNOSIS, levels = c('HC','MCI','AD'))

ggplot(MTecc, aes(x = POSITION, y = MT, group = DIAGNOSIS, colour = DIAGNOSIS,
                   shape = DIAGNOSIS)) +
  geom_point(size = 3, position = position_dodge(width = .4)) +
  geom_errorbar(aes(ymin=MT-ci, ymax=MT+ci), 
                width=.4, position = position_dodge(width = .4)) + 
  geom_line(aes(group = DIAGNOSIS), size = 0.7, position = position_dodge(width = .4)) +
  labs(x = 'Eccentricity (°)', y = 'Movement time (ms)') +
  scale_color_manual(values = c('black','grey30','grey60')) +
  facet_wrap(~VIEW) + theme_classic() + ylim(500,1000) +
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12)
  )

ggsave('RADMTecc_plot.png', plot = last_plot(),  device = NULL, dpi = 300, 
       width = 7.5, height = 5, path = anaPath)

## ANOVA ## 
MTanova <- MT_means[MT_means$PPT != 212 & MT_means$PPT != 310 ,]

# FULL ANOVA ON MOVEMENT TIME DATA
MT_ANOVA <- ezANOVA(
  data = MTanova
  , dv = .(MT)
  , wid = .(PPT)
  , within = .(VIEW, DOM)
  , between = .(DIAGNOSIS)
  , between_covariates = .(SITE)
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
# removing ppt with fewer trials (data invalid)
MTECCanova <- MT_medians[MT_medians$PPT != 212 & 
                           MT_medians$PPT != 302 & 
                           MT_medians$PPT != 310 &
                           MT_medians$PPT != 315 &
                           MT_medians$PPT != 403 &
                           MT_medians$PPT != 407,]
MTECCanova$ECC <- factor(MTECCanova$ECC)

MTECCanova$POSITION <- factor(MTECCanova$POSITION)
MTECC_ANOVA <- ezANOVA(
  data = MTECCanova
  , dv = .(MT)
  , wid = .(PPT)
  , within = .(VIEW, DOM, ECC)
  , between = .(DIAGNOSIS)
  , between_covariates = .(SITE)
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
RT_medians <- aggregate(RT ~ PPT * VIEW * SIDE * DOM * POSITION * ECC * SITE * GRP * DIAGNOSIS, 
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

ggsave('RADRTmean_plot.png', plot = last_plot(),  device = NULL, dpi = 300, 
       width = 4, height = 5, path = anaPath)

# average across sides
RTav$DIAGNOSIS <- factor(RTav$DIAGNOSIS, levels = c('HC','MCI','AD'))
ggplot(RTav, aes(x = VIEW, y = RT, colour = SITE, group = PPT)) +
  geom_point(shape = 16, size = 2, position = position_dodge(width = .3)) +
  geom_line(aes(group = PPT), size = 0.5, alpha = .5, 
            position = position_dodge(width = .3)) +
  scale_color_manual(values = c('dodgerblue3','grey50')) +
  stat_summary(aes(y = RT, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1, size = 4, group = 1) +
  facet_wrap(~DIAGNOSIS) + 
  labs(x = '', y = 'Reaction time (ms)', element_text(size = 12)) +
  theme_classic() + theme(legend.position = 'bottom', 
                          axis.text = element_text(size = 10),
                          axis.title = element_text(size = 12),
                          strip.text = element_text(size = 12)
  ) 

ggsave('RADRTav_plot.png', plot = last_plot(),  device = NULL, dpi = 300, 
       width = 6, height = 4, path = anaPath)

# by eccentricity
RT_medians$POSITION <- factor(RT_medians$POSITION)

# summary data
RTecc <- summarySE(RT_medians, measurevar = 'RT', 
                   groupvar = c('DIAGNOSIS','POSITION','VIEW','SIDE'), na.rm = TRUE)
RTecc$DIAGNOSIS <- factor(RTecc$DIAGNOSIS, levels = c('HC','MCI','AD'))

ggplot(RTecc, aes(x = POSITION, y = RT, group = DIAGNOSIS, colour = DIAGNOSIS,
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


ggsave('RADRTecc_plot.png', plot = last_plot(),  device = NULL, dpi = 300, 
       width = 7.5, height = 5, path = anaPath)

## ANOVA ##
# FULL ANOVA ON REACTION TIME DATA BY ECCENTRICITY
RTECCanova <- RT_medians[RT_medians$PPT != 212 & 
                           RT_medians$PPT != 310 &
                           RT_medians$PPT != 315 &
                           RT_medians$PPT != 302 &
                           RT_medians$PPT != 403 &
                           RT_medians$PPT != 407,]
RTECCanova$ECC <- factor(RTECCanova$ECC)

RTECC_ANOVA <- ezANOVA(
  data = RTECCanova
  , dv = .(RT)
  , wid = .(PPT)
  , within = .(VIEW, DOM, ECC)
  , between = .(DIAGNOSIS)
  , between_covariates = .(SITE)
  , type = 3,
  return_aov = TRUE,
  detailed = TRUE
)

RTECC_ANOVA$ANOVA
RTECC_ANOVA$`Mauchly's Test for Sphericity`
RTECC_ANOVA$`Sphericity Corrections`
aovRTECC <- aovEffectSize(ezObj = RTECC_ANOVA, effectSize = "pes")
aovDispTable(aovRTECC)

#pair-wise t-test
RTttest <- pairwise.t.test(RT_medians$RT, RT_medians$DIAGNOSIS, p.adj = 'bonf')
print(RTttest)

##### PEAK SPEED #####
PS_medians <- aggregate(PS ~ PPT * VIEW * SIDE * DOM * POSITION * ECC * SITE * GRP * DIAGNOSIS, 
                        median, data = res)
# remove 310, not using them and 315 - noisy
PS_medians <- PS_medians[PS_medians$PPT != 310 & PS_medians$PPT != 315 ,]
# remove outlier in P407
PS_medians <- PS_medians[!(PS_medians$PPT == '407' & PS_medians$PS > 4000) ,]
PS_means <- aggregate(PS ~ PPT * VIEW * SIDE * DOM * SITE * GRP * DIAGNOSIS,
                      mean, data = PS_medians)
# average across side
PSav <- aggregate(PS ~ PPT * VIEW * SITE * GRP * DIAGNOSIS,
                  mean, data = PS_medians)

## plotting :)
# both sides
PSsummary <- summarySE(PS_means, measurevar = 'PS', 
                       groupvar = c('DIAGNOSIS','VIEW','DOM'), na.rm = TRUE)
PSsummary$DIAGNOSIS <- factor(PSsummary$DIAGNOSIS, levels = c('HC','MCI','AD'))
PSsummary$DOM <- factor(PSsummary$DOM, labels = c('Non-dominant','Dominant'))

ggplot(PSsummary, aes(x = VIEW, y = PS, colour = DIAGNOSIS, group = DIAGNOSIS,
                      shape = DIAGNOSIS)) +
  geom_point(size = 3, position = position_dodge(width = .3)) +
  geom_line(aes(group = DIAGNOSIS), size = 0.5, alpha = .5,
            position = position_dodge(width = .3)) +
  geom_errorbar(aes(ymin=PS-ci, ymax=PS+ci), 
                width=.4, position = position_dodge(width = .3)) +
  scale_color_manual(values = c('black','grey30','grey60')) +
  facet_grid(~DOM) + #ylim(300,900) +
  labs(x = 'Side', y = 'Peak velocity (mm/s2)', element_text(size = 12)) +
  theme_classic() + theme(legend.position = 'bottom', 
                          legend.title = element_blank(),
                          axis.text = element_text(size = 10),
                          axis.title = element_text(size = 12),
                          strip.text = element_text(size = 12)
  ) 

ggsave('RADPSmean_plot.png', plot = last_plot(),  device = NULL, dpi = 300, 
       width = 4, height = 5, path = anaPath)

# average across sides
PSav$DIAGNOSIS <- factor(PSav$DIAGNOSIS, levels = c('HC','MCI','AD'))
ggplot(PSav, aes(x = VIEW, y = PS, colour = SITE, group = PPT)) +
  geom_point(shape = 16, size = 2, position = position_dodge(width = .3)) +
  geom_line(aes(group = PPT), size = 0.5, alpha = .5, 
            position = position_dodge(width = .3)) +
  scale_color_manual(values = c('dodgerblue3','grey50')) +
  stat_summary(aes(y = PS, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1, size = 4, group = 1) +
  facet_wrap(~DIAGNOSIS) + 
  labs(x = '', y = 'Peak velocity (mm/s2)', element_text(size = 12)) +
  theme_classic() + theme(legend.position = 'bottom', 
                          axis.text = element_text(size = 10),
                          axis.title = element_text(size = 12),
                          strip.text = element_text(size = 12)
  ) 

ggsave('RADPSav_plot.png', plot = last_plot(),  device = NULL, dpi = 300, 
       width = 6, height = 4, path = anaPath)

# by eccentricity
PS_medians$POSITION <- factor(PS_medians$POSITION)

# summary data
PSecc <- summarySE(PS_medians, measurevar = 'PS', 
                   groupvar = c('DIAGNOSIS','POSITION','VIEW','SIDE'), na.rm = TRUE)
PSecc$DIAGNOSIS <- factor(PSecc$DIAGNOSIS, levels = c('HC','MCI','AD'))

ggplot(PSecc, aes(x = POSITION, y = PS, group = DIAGNOSIS, colour = DIAGNOSIS,
                  shape = DIAGNOSIS)) +
  geom_point(size = 3, position = position_dodge(width = .4)) +
  geom_errorbar(aes(ymin=PS-ci, ymax=PS+ci), 
                width=.4, position = position_dodge(width = .4)) + 
  geom_line(aes(group = DIAGNOSIS), size = 0.7, position = position_dodge(width = .4)) +
  labs(x = 'Eccentricity (mm)', y = 'Peak velocity (mm/s2)') +
  scale_color_manual(values = c('black','grey30','grey60')) +
  facet_grid(~VIEW) + theme_classic() + 
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12)
  )


ggsave('RADPSecc_plot.png', plot = last_plot(),  device = NULL, dpi = 300, 
       width = 7.5, height = 5, path = anaPath)

## ANOVA ##
# full anova on PS by ecc
PSECCanova <- PS_medians[PS_medians$PPT != 212 & 
                           PS_medians$PPT != 403 &
                           PS_medians$PPT != 407 ,]
PSECCanova$ECC <- factor(PSECCanova$ECC)

PSECC_ANOVA <- ezANOVA(
  data = PSECCanova
  , dv = .(PS)
  , wid = .(PPT)
  , within = .(VIEW, DOM, ECC)
  , between = .(DIAGNOSIS)
  , between_covariates = .(SITE)
  , type = 3,
  return_aov = TRUE,
  detailed = TRUE
)

PSECC_ANOVA$ANOVA
PSECC_ANOVA$`Mauchly's Test for Sphericity`
PSECC_ANOVA$`Sphericity Corrections`
aovPSECC <- aovEffectSize(ezObj = PSECC_ANOVA, effectSize = "pes")
aovDispTable(aovPSECC)

#pair-wise t-test
PSttest <- pairwise.t.test(PS_medians$PS, PS_medians$DIAGNOSIS, p.adj = 'bonf')
print(PSttest)

##### TIME TO PV #####
TPS_medians <- aggregate(TPS ~ PPT * VIEW * SIDE * DOM * POSITION * ECC * SITE * GRP * DIAGNOSIS, 
                         median, data = res)
# remove 310, not using them and 315 - noisy
TPS_medians <- TPS_medians[TPS_medians$PPT != 310 & TPS_medians$PPT != 315 ,]
TPS_means <- aggregate(TPS ~ PPT * VIEW * SIDE * DOM * SITE * GRP * DIAGNOSIS,
                        mean, data = TPS_medians)
# average across side
TPSav <- aggregate(TPS ~ PPT * VIEW * SITE * GRP * DIAGNOSIS,
                    mean, data = TPS_medians)

## plotting :)
# both sides
TPSsummary <- summarySE(TPS_means, measurevar = 'TPS', 
                         groupvar = c('DIAGNOSIS','VIEW','DOM'), na.rm = TRUE)
TPSsummary$DIAGNOSIS <- factor(TPSsummary$DIAGNOSIS, levels = c('HC','MCI','AD'))
TPSsummary$DOM <- factor(TPSsummary$DOM, labels = c('Non-dominant','Dominant'))

ggplot(TPSsummary, aes(x = VIEW, y = TPS, colour = DIAGNOSIS, group = DIAGNOSIS,
                        shape = DIAGNOSIS)) +
  geom_point(size = 3, position = position_dodge(width = .3)) +
  geom_line(aes(group = DIAGNOSIS), size = 0.5, alpha = .5,
            position = position_dodge(width = .3)) +
  geom_errorbar(aes(ymin=TPS-ci, ymax=TPS+ci), 
                width=.4, position = position_dodge(width = .3)) +
  scale_color_manual(values = c('black','grey30','grey60')) +
  facet_grid(~DOM) + 
  labs(x = 'Side', y = 'Time to peak velocity (ms)', element_text(size = 12)) +
  theme_classic() + theme(legend.position = 'bottom', 
                          legend.title = element_blank(),
                          axis.text = element_text(size = 10),
                          axis.title = element_text(size = 12),
                          strip.text = element_text(size = 12)
  ) 

ggsave('RADTPSmean_plot.png', plot = last_plot(),  device = NULL, dpi = 300, 
       width = 4, height = 5, path = anaPath)

# average across sides
TPSav$DIAGNOSIS <- factor(TPSav$DIAGNOSIS, levels = c('HC','MCI','AD'))
ggplot(TPSav, aes(x = VIEW, y = TPS, colour = SITE, group = PPT)) +
  geom_point(shape = 16, size = 2, position = position_dodge(width = .3)) +
  geom_line(aes(group = PPT), size = 0.5, alpha = .5, 
            position = position_dodge(width = .3)) +
  scale_color_manual(values = c('dodgerblue3','grey50')) +
  stat_summary(aes(y = TPS, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1, size = 4, group = 1) +
  facet_wrap(~DIAGNOSIS) + 
  labs(x = '', y = 'Time to peak velocity (ms)', element_text(size = 12)) +
  theme_classic() + theme(legend.position = 'bottom', 
                          axis.text = element_text(size = 10),
                          axis.title = element_text(size = 12),
                          strip.text = element_text(size = 12)
  ) 

ggsave('RADTPSav_plot.png', plot = last_plot(),  device = NULL, dpi = 300, 
       width = 6, height = 4, path = anaPath)

# by eccentricity
TPS_medians$POSITION <- factor(TPS_medians$POSITION)

# summary data
TPSecc <- summarySE(TPS_medians, measurevar = 'TPS', 
                     groupvar = c('DIAGNOSIS','POSITION','VIEW','SIDE'), na.rm = TRUE)
TPSecc$DIAGNOSIS <- factor(TPSecc$DIAGNOSIS, levels = c('HC','MCI','AD'))

ggplot(TPSecc, aes(x = POSITION, y = TPS, group = DIAGNOSIS, colour = DIAGNOSIS,
                    shape = DIAGNOSIS)) +
  geom_point(size = 3, position = position_dodge(width = .4)) +
  geom_errorbar(aes(ymin=TPS-ci, ymax=TPS+ci), 
                width=.4, position = position_dodge(width = .4)) + 
  geom_line(aes(group = DIAGNOSIS), size = 0.7, position = position_dodge(width = .4)) +
  labs(x = 'Eccentricity (mm)', y = 'Time to peak velocity (ms)') +
  scale_color_manual(values = c('black','grey30','grey60')) +
  facet_grid(~VIEW) + theme_classic() + ylim(150,280) +
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12)
  )


ggsave('RADTPSecc_plot.png', plot = last_plot(),  device = NULL, dpi = 300, 
       width = 7.5, height = 5, path = anaPath)

## ANOVA ##
# full anova on TPS by ecc
TPSECCanova <- TPS_medians[TPS_medians$PPT != 212 & 
                             TPS_medians$PPT != 302 &
                             TPS_medians$PPT != 403 &
                             TPS_medians$PPT != 407 ,]
TPSECCanova$ECC <- factor(TPSECCanova$ECC)

TPSECC_ANOVA <- ezANOVA(
  data = TPSECCanova
  , dv = .(TPS)
  , wid = .(PPT)
  , within = .(VIEW, DOM, ECC)
  , between = .(DIAGNOSIS)
  , between_covariates = .(SITE)
  , type = 3,
  return_aov = TRUE,
  detailed = TRUE
)

TPSECC_ANOVA$ANOVA
TPSECC_ANOVA$`Mauchly's Test for Sphericity`
TPSECC_ANOVA$`Sphericity Corrections`
aovTPSECC <- aovEffectSize(ezObj = TPSECC_ANOVA, effectSize = "pes")
aovDispTable(aovTPSECC)

#pair-wise t-test
TPSttest <- pairwise.t.test(TPS_medians$TPS, TPS_medians$DIAGNOSIS, p.adj = 'bonf')
print(TPSttest)

##### TIME AFTER PV #####
res$TAPS <- res$MT - res$TPS

TAPS_medians <- aggregate(TAPS ~ PPT * VIEW * SIDE * DOM * POSITION * ECC * SITE * GRP * DIAGNOSIS, 
                          median, data = res)
# remove 310, not using them and 315 - noisy
TAPS_medians <- TAPS_medians[TAPS_medians$PPT != 310 & TAPS_medians$PPT != 315 ,]
TAPS_means <- aggregate(TAPS ~ PPT * VIEW * SIDE * DOM * SITE * GRP * DIAGNOSIS,
                        mean, data = TAPS_medians)
# average across side
TAPSav <- aggregate(TAPS ~ PPT * VIEW * SITE * GRP * DIAGNOSIS,
                    mean, data = TAPS_medians)

## plotting :)
# both sides
TAPSsummary <- summarySE(TAPS_means, measurevar = 'TAPS', 
                         groupvar = c('DIAGNOSIS','VIEW','DOM'), na.rm = TRUE)
TAPSsummary$DIAGNOSIS <- factor(TAPSsummary$DIAGNOSIS, levels = c('HC','MCI','AD'))
TAPSsummary$DOM <- factor(TAPSsummary$DOM, labels = c('Non-dominant','Dominant'))

ggplot(TAPSsummary, aes(x = VIEW, y = TAPS, colour = DIAGNOSIS, group = DIAGNOSIS,
                        shape = DIAGNOSIS)) +
  geom_point(size = 3, position = position_dodge(width = .3)) +
  geom_line(aes(group = DIAGNOSIS), size = 0.5, alpha = .5,
            position = position_dodge(width = .3)) +
  geom_errorbar(aes(ymin=TAPS-ci, ymax=TAPS+ci), 
                width=.4, position = position_dodge(width = .3)) +
  scale_color_manual(values = c('black','grey30','grey60')) +
  facet_grid(~DOM) + ylim(300,600) +
  labs(x = 'Side', y = 'Time after peak velocity (ms)', element_text(size = 12)) +
  theme_classic() + theme(legend.position = 'bottom', 
                          legend.title = element_blank(),
                          axis.text = element_text(size = 10),
                          axis.title = element_text(size = 12),
                          strip.text = element_text(size = 12)
  ) 

ggsave('RADTAPSmean_plot.png', plot = last_plot(),  device = NULL, dpi = 300, 
       width = 4, height = 5, path = anaPath)

# average across sides
TAPSav$DIAGNOSIS <- factor(TAPSav$DIAGNOSIS, levels = c('HC','MCI','AD'))
ggplot(TAPSav, aes(x = VIEW, y = TAPS, colour = SITE, group = PPT)) +
  geom_point(shape = 16, size = 2, position = position_dodge(width = .3)) +
  geom_line(aes(group = PPT), size = 0.5, alpha = .5, 
            position = position_dodge(width = .3)) +
  scale_color_manual(values = c('dodgerblue3','grey50')) +
  stat_summary(aes(y = TAPS, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1, size = 4, group = 1) +
  facet_wrap(~DIAGNOSIS) + #ylim(0.5,0.8) +
  labs(x = '', y = 'Time after peak velocity (ms)', element_text(size = 12)) +
  theme_classic() + theme(legend.position = 'bottom', 
                          axis.text = element_text(size = 10),
                          axis.title = element_text(size = 12),
                          strip.text = element_text(size = 12)
  ) 

ggsave('RADTAPSav_plot.png', plot = last_plot(),  device = NULL, dpi = 300, 
       width = 6, height = 4, path = anaPath)

# by eccentricity
TAPS_medians$POSITION <- factor(TAPS_medians$POSITION)

# summary data
TAPSecc <- summarySE(TAPS_medians, measurevar = 'TAPS', 
                     groupvar = c('DIAGNOSIS','POSITION','VIEW','SIDE'), na.rm = TRUE)
TAPSecc$DIAGNOSIS <- factor(TAPSecc$DIAGNOSIS, levels = c('HC','MCI','AD'))

ggplot(TAPSecc, aes(x = POSITION, y = TAPS, group = DIAGNOSIS, colour = DIAGNOSIS,
                    shape = DIAGNOSIS)) +
  geom_point(size = 3, position = position_dodge(width = .4)) +
  geom_errorbar(aes(ymin=TAPS-ci, ymax=TAPS+ci), 
                width=.4, position = position_dodge(width = .4)) + 
  geom_line(aes(group = DIAGNOSIS), size = 0.7, position = position_dodge(width = .4)) +
  labs(x = 'Eccentricity (mm)', y = 'Time after peak velocity (ms)') +
  scale_color_manual(values = c('black','grey30','grey60')) +
  facet_grid(~VIEW) + theme_classic() + #ylim(0.6,0.8) +
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12)
  )


ggsave('RADTAPSecc_plot.png', plot = last_plot(),  device = NULL, dpi = 300, 
       width = 7.5, height = 5, path = anaPath)

## ANOVA ##
# full anova on TAPS by ecc
TAPSECCanova <- TAPS_medians[TAPS_medians$PPT != 212 & 
                               TAPS_medians$PPT != 302 &
                               TAPS_medians$PPT != 403 &
                               TAPS_medians$PPT != 407,]
TAPSECCanova$ECC <- factor(TAPSECCanova$ECC)

TAPSECC_ANOVA <- ezANOVA(
  data = TAPSECCanova
  , dv = .(TAPS)
  , wid = .(PPT)
  , within = .(VIEW, DOM, ECC)
  , between = .(DIAGNOSIS)
  , between_covariates = .(SITE)
  , type = 3,
  return_aov = TRUE,
  detailed = TRUE
)

TAPSECC_ANOVA$ANOVA
TAPSECC_ANOVA$`Mauchly's Test for Sphericity`
TAPSECC_ANOVA$`Sphericity Corrections`
aovTAPSECC <- aovEffectSize(ezObj = TAPSECC_ANOVA, effectSize = "pes")
aovDispTable(aovTAPSECC)

#pair-wise t-test
TAPSttest <- pairwise.t.test(TAPS_medians$TAPS, TAPS_medians$DIAGNOSIS, p.adj = 'bonf')
print(TAPSttest)

##### NORM TIME AFTER PV #####
res$NTAPS <- (res$MT - res$TPS)/res$MT

NTAPS_medians <- aggregate(NTAPS ~ PPT * VIEW * SIDE * DOM * POSITION * ECC * SITE * GRP * DIAGNOSIS, 
                        median, data = res)
# remove 310, not using them and 315 - noisy
NTAPS_medians <- NTAPS_medians[NTAPS_medians$PPT != 310 & NTAPS_medians$PPT != 315 ,]
NTAPS_means <- aggregate(NTAPS ~ PPT * VIEW * SIDE * DOM * SITE * GRP * DIAGNOSIS,
                      mean, data = NTAPS_medians)
# average across side
NTAPSav <- aggregate(NTAPS ~ PPT * VIEW * SITE * GRP * DIAGNOSIS,
                  mean, data = NTAPS_medians)

## plotting :)
# both sides
NTAPSsummary <- summarySE(NTAPS_means, measurevar = 'NTAPS', 
                       groupvar = c('DIAGNOSIS','VIEW','DOM'), na.rm = TRUE)
NTAPSsummary$DIAGNOSIS <- factor(NTAPSsummary$DIAGNOSIS, levels = c('HC','MCI','AD'))
NTAPSsummary$DOM <- factor(NTAPSsummary$DOM, labels = c('Non-dominant','Dominant'))

ggplot(NTAPSsummary, aes(x = VIEW, y = NTAPS, colour = DIAGNOSIS, group = DIAGNOSIS,
                      shape = DIAGNOSIS)) +
  geom_point(size = 3, position = position_dodge(width = .3)) +
  geom_line(aes(group = DIAGNOSIS), size = 0.5, alpha = .5,
            position = position_dodge(width = .3)) +
  geom_errorbar(aes(ymin=NTAPS-ci, ymax=NTAPS+ci), 
                width=.4, position = position_dodge(width = .3)) +
  scale_color_manual(values = c('black','grey30','grey60')) +
  facet_grid(~DOM) + ylim(0.6,0.8) +
  labs(x = 'Side', y = 'Normalised time after peak velocity', element_text(size = 12)) +
  theme_classic() + theme(legend.position = 'bottom', 
                          legend.title = element_blank(),
                          axis.text = element_text(size = 10),
                          axis.title = element_text(size = 12),
                          strip.text = element_text(size = 12)
  ) 

ggsave('RADnTAPSmean_plot.png', plot = last_plot(),  device = NULL, dpi = 300, 
       width = 4, height = 5, path = anaPath)

# average across sides
NTAPSav$DIAGNOSIS <- factor(NTAPSav$DIAGNOSIS, levels = c('HC','MCI','AD'))
ggplot(NTAPSav, aes(x = VIEW, y = NTAPS, colour = SITE, group = PPT)) +
  geom_point(shape = 16, size = 2, position = position_dodge(width = .3)) +
  geom_line(aes(group = PPT), size = 0.5, alpha = .5, 
            position = position_dodge(width = .3)) +
  scale_color_manual(values = c('dodgerblue3','grey50')) +
  stat_summary(aes(y = NTAPS, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1, size = 4, group = 1) +
  facet_wrap(~DIAGNOSIS) + ylim(0.5,0.8) +
  labs(x = '', y = 'Normalised time after peak velocity', element_text(size = 12)) +
  theme_classic() + theme(legend.position = 'bottom', 
                          axis.text = element_text(size = 10),
                          axis.title = element_text(size = 12),
                          strip.text = element_text(size = 12)
  ) 

ggsave('RADnTAPSav_plot.png', plot = last_plot(),  device = NULL, dpi = 300, 
       width = 6, height = 4, path = anaPath)

# by eccentricity
NTAPS_medians$POSITION <- factor(NTAPS_medians$POSITION)

# summary data
NTAPSecc <- summarySE(NTAPS_medians, measurevar = 'NTAPS', 
                   groupvar = c('DIAGNOSIS','POSITION','VIEW','SIDE'), na.rm = TRUE)
NTAPSecc$DIAGNOSIS <- factor(NTAPSecc$DIAGNOSIS, levels = c('HC','MCI','AD'))

ggplot(NTAPSecc, aes(x = POSITION, y = NTAPS, group = DIAGNOSIS, colour = DIAGNOSIS,
                  shape = DIAGNOSIS)) +
  geom_point(size = 3, position = position_dodge(width = .4)) +
  geom_errorbar(aes(ymin=NTAPS-ci, ymax=NTAPS+ci), 
                width=.4, position = position_dodge(width = .4)) + 
  geom_line(aes(group = DIAGNOSIS), size = 0.7, position = position_dodge(width = .4)) +
  labs(x = 'Eccentricity (mm)', y = 'Normalised time after peak velocity') +
  scale_color_manual(values = c('black','grey30','grey60')) +
  facet_grid(~VIEW) + theme_classic() + ylim(0.6,0.8) +
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12)
  )


ggsave('RADnTAPSecc_plot.png', plot = last_plot(),  device = NULL, dpi = 300, 
       width = 7.5, height = 5, path = anaPath)

## ANOVA ##
# full anova on TAPS by ecc
NTAPSECCanova <- NTAPS_medians[NTAPS_medians$PPT != 212 & 
                                 NTAPS_medians$PPT != 302 &
                                 NTAPS_medians$PPT != 403 &
                                 NTAPS_medians$PPT != 407,]
NTAPSECCanova$ECC <- factor(NTAPSECCanova$ECC)

NTAPSECC_ANOVA <- ezANOVA(
  data = NTAPSECCanova
  , dv = .(NTAPS)
  , wid = .(PPT)
  , within = .(VIEW, DOM, ECC)
  , between = .(DIAGNOSIS)
  , between_covariates = .(SITE)
  , type = 3,
  return_aov = TRUE,
  detailed = TRUE
)

NTAPSECC_ANOVA$ANOVA
NTAPSECC_ANOVA$`Mauchly's Test for Sphericity`
NTAPSECC_ANOVA$`Sphericity Corrections`
aovNTAPSECC <- aovEffectSize(ezObj = NTAPSECC_ANOVA, effectSize = "pes")
aovDispTable(aovNTAPSECC)

#pair-wise t-test
NTAPSttest <- pairwise.t.test(NTAPS_medians$NTAPS, NTAPS_medians$DIAGNOSIS, p.adj = 'bonf')
print(NTAPSttest)

##### CORRELATE ACE #####
setwd(dataPath)
patient_demos <- read.csv('patient_demographics.csv')
PMIdata <- read.csv('radialPMI-filtered.csv')
#extracting ACE data into seperate data-frame
ACEscores <- patient_demos[ ,c(1, 10:15)]
names(ACEscores)[1] <- 'PPT'
setwd(anaPath)
# merging ACE with PMI
PMIACE <- merge(PMIdata, ACEscores, by = 'PPT')
PMIACE$ACEall <- as.numeric(as.character(PMIACE$ACEall))
PMIACE$ACEvisuospatial <- as.numeric(as.character(PMIACE$ACEvisuospatial))

## aaaand correlate :)
ggscatter(PMIACE, x = 'ACEall', y = 'PMI', add = 'reg.line', conf.int = TRUE,
          cor.coef = TRUE, size = 1, cor.coef.size = 3, cor.method = 'spearman') +
  ylab('Radial PMI (deg)') + xlab('ACE score (%)') +
  facet_wrap(~DOM) +
  theme(text = element_text(size = 10))

# average PMI
PMIACE <- aggregate(PMI ~ PPT+DIAGNOSIS+ACEall, mean, data = PMIACE)
ggscatter(PMIACE, x = 'ACEall', y = 'PMI', 
          add = 'reg.line', conf.int = FALSE, add.params = list(color = "black"),
          cor.coef = TRUE, size = 1.5, cor.coef.size = 3, cor.method = 'spearman') +
  ylab('Radial PMI (deg)') + xlab('ACE score (%)') +
  theme(text = element_text(size = 10))

ggsave('RADPMI-ACEcorr_plot.png', plot = last_plot(),  device = NULL, dpi = 300, 
       width = 5, height = 5, path = anaPath)

