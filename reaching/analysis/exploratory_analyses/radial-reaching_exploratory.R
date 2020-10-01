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
ggplot(ANGmeans, aes(x = VIEW, y = ANG_ERR, group = PPT, colour = SITE)) +
  geom_hline(yintercept = 0) +
  geom_point(shape = 16, size = 2, position = position_dodge(width = .3)) +
  stat_summary(aes(y = ANG_ERR, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1, size = 4, group = 1) +
  geom_line(aes(group = PPT), alpha = .4, size = 0.7, position = position_dodge(width = .3)) +
  scale_color_manual(values = c('dodgerblue3','grey50')) + 
  facet_grid(~DIAGNOSIS) +
  labs(x = 'Side', y = 'Angular error (°)') + 
  theme_classic() + theme(legend.position = 'bottom',
                          axis.text = element_text(size = 10),
                          axis.title = element_text(size = 12),
                          strip.text = element_text(size = 12)
                          )

ggsave('RAD_ANGERRmeans.png', plot = last_plot(), device = NULL, dpi = 300, 
       width = 5, height = 7, path = anaPath)

## eccentricity
# plotting data frame
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
  labs(x = 'Eccentricity (mm)', y = 'Angular error (°)') + 
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
ggplot(AMPmeans, aes(x = VIEW, y = AMP_ERR, group = PPT, colour = SITE)) +
  geom_hline(yintercept = 0) +
  geom_point(shape = 16, size = 2, position = position_dodge(width = .3)) +
  stat_summary(aes(y = AMP_ERR, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1, size = 4, group = 1) +
  geom_line(aes(group = PPT), alpha = .4, size = 0.7, position = position_dodge(width = .3)) +
  scale_color_manual(values = c('dodgerblue3','grey50')) +
  facet_grid(~DIAGNOSIS) +
  labs(x = 'Side', y = 'Amplutude error (mm)') + 
  theme_classic() + theme(legend.position = 'bottom',
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
  labs(x = 'Eccentricity (mm)', y = 'Amplitude error (mm)') + 
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
MT_ECC <- aggregate(MT ~ PPT * VIEW * ECC * SITE * DIAGNOSIS, 
                        median, data = MT_medians)
MT_means <- aggregate(MT ~ PPT * VIEW * SITE * DIAGNOSIS,
                      mean, data = MT_ECC)

## plotting :)
# means for all PP
MT_means$DIAGNOSIS <- factor(MT_means$DIAGNOSIS, levels = c('HC','MCI','AD'))

ggplot(MT_means, aes(x = VIEW, y = MT, colour = SITE, group = PPT)) +
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

ggsave('RAD_MTmeans.png', plot = last_plot(),  device = NULL, dpi = 300, 
       width = 6, height = 4, path = anaPath)

## by eccentricity
# summary data
MTecc <- summarySE(MT_ECC, measurevar = 'MT', 
                   groupvar = c('DIAGNOSIS','ECC','VIEW'), na.rm = TRUE)
MTecc$DIAGNOSIS <- factor(MTecc$DIAGNOSIS, levels = c('HC','MCI','AD'))

ggplot(MTecc, aes(x = ECC, y = MT, group = DIAGNOSIS, colour = DIAGNOSIS,
                   shape = DIAGNOSIS)) +
  geom_point(size = 3, position = position_dodge(width = .4)) +
  geom_errorbar(aes(ymin=MT-ci, ymax=MT+ci), 
                width=.4, position = position_dodge(width = .4)) + 
  geom_line(aes(group = DIAGNOSIS), size = 0.7, position = position_dodge(width = .4)) +
  labs(x = 'Eccentricity (mm)', y = 'Movement time (ms)') +
  scale_color_manual(values = c('black','grey30','grey60')) +
  facet_wrap(~VIEW) + theme_classic() + ylim(500,1000) +
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12)
  )

ggsave('RAD_MT_ECC.png', plot = last_plot(),  device = NULL, dpi = 300, 
       width = 7.5, height = 5, path = anaPath)

## ANOVA ## 
# FULL ANOVA ON MOVEMENT TIME DATA
MT_ANOVA <- ezANOVA(
  data = MT_ECC
  , dv = .(MT)
  , wid = .(PPT)
  , within = .(VIEW, ECC)
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

#pair-wise t-test
MTttest <- pairwise.t.test(MT_ECC$MT, MT_ECC$DIAGNOSIS, p.adj = 'bonf')
print(MTttest)

###### REACTION TIME ######
RT_medians <- aggregate(RT ~ PPT * VIEW * SIDE * DOM * ECC * SITE * DIAGNOSIS, 
                        median, data = res)
RT_ECC <- aggregate(RT ~ PPT * VIEW * ECC * SITE * DIAGNOSIS, 
                    median, data = RT_medians)
RT_means <- aggregate(RT ~ PPT * VIEW * SITE * DIAGNOSIS,
                      mean, data = RT_ECC)

## plotting 
# mean - all PP
RT_means$DIAGNOSIS <- factor(RT_means$DIAGNOSIS, levels = c('HC','MCI','AD'))

ggplot(RT_means, aes(x = VIEW, y = RT, colour = SITE, group = PPT)) +
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

ggsave('RAD_RTmean.png', plot = last_plot(),  device = NULL, dpi = 300, 
       width = 6, height = 4, path = anaPath)

## by eccentricity
# summary data
RTecc <- summarySE(RT_ECC, measurevar = 'RT', 
                   groupvar = c('DIAGNOSIS','ECC','VIEW'), na.rm = TRUE)
RTecc$DIAGNOSIS <- factor(RTecc$DIAGNOSIS, levels = c('HC','MCI','AD'))

ggplot(RTecc, aes(x = ECC, y = RT, group = DIAGNOSIS, colour = DIAGNOSIS,
                 shape = DIAGNOSIS)) +
  geom_point(size = 3, position = position_dodge(width = .4)) +
  geom_errorbar(aes(ymin=RT-ci, ymax=RT+ci), 
                width=.4, position = position_dodge(width = .4)) + 
  geom_line(aes(group = DIAGNOSIS), size = 0.7, position = position_dodge(width = .4)) +
  labs(x = 'Eccentricity (mm)', y = 'Reaction time (ms)') +
  scale_color_manual(values = c('black','grey30','grey60')) +
  facet_grid(~VIEW) + theme_classic() +
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12)
  )


ggsave('RAD_RT_ECC.png', plot = last_plot(),  device = NULL, dpi = 300, 
       width = 7.5, height = 5, path = anaPath)

## ANOVA ##
RTECC_ANOVA <- ezANOVA(
  data = RT_ECC
  , dv = .(RT)
  , wid = .(PPT)
  , within = .(VIEW, ECC)
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
RTttest <- pairwise.t.test(RT_ECC$RT, RT_ECC$DIAGNOSIS, p.adj = 'bonf')
print(RTttest)

##### PEAK SPEED #####
PS_medians <- aggregate(PS ~ PPT * VIEW * SIDE * DOM * POSITION * ECC * SITE * DIAGNOSIS, 
                        median, data = res)
PS_ECC <- aggregate(PS ~ PPT * VIEW * ECC * SITE * DIAGNOSIS, 
                        median, data = PS_medians)
PS_means <- aggregate(PS ~ PPT * VIEW * SITE * DIAGNOSIS,
                      mean, data = PS_ECC)

## plotting :)
# mean all PP
PS_means$DIAGNOSIS <- factor(PS_means$DIAGNOSIS, levels = c('HC','MCI','AD'))

ggplot(PS_means, aes(x = VIEW, y = PS, colour = SITE, group = PPT)) +
  geom_point(shape = 16, size = 2, position = position_dodge(width = .3)) +
  geom_line(aes(group = PPT), size = 0.5, alpha = .5, 
            position = position_dodge(width = .3)) +
  scale_color_manual(values = c('dodgerblue3','grey50')) +
  stat_summary(aes(y = PS, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1, size = 4, group = 1) +
  facet_wrap(~DIAGNOSIS) + 
  labs(x = '', y = 'Peak speed (ms)', element_text(size = 12)) +
  theme_classic() + theme(legend.position = 'bottom', 
                          axis.text = element_text(size = 10),
                          axis.title = element_text(size = 12),
                          strip.text = element_text(size = 12)
  ) 

ggsave('RAD_PSmean.png', plot = last_plot(),  device = NULL, dpi = 300, 
       width = 6, height = 4, path = anaPath)

# by eccentricity
# summary data
PSecc <- summarySE(PS_medians, measurevar = 'PS', 
                   groupvar = c('DIAGNOSIS','ECC','VIEW'), na.rm = TRUE)
PSecc$DIAGNOSIS <- factor(PSecc$DIAGNOSIS, levels = c('HC','MCI','AD'))

ggplot(PSecc, aes(x = ECC, y = PS, group = DIAGNOSIS, colour = DIAGNOSIS,
                  shape = DIAGNOSIS)) +
  geom_point(size = 3, position = position_dodge(width = .4)) +
  geom_errorbar(aes(ymin=PS-ci, ymax=PS+ci), 
                width=.4, position = position_dodge(width = .4)) + 
  geom_line(aes(group = DIAGNOSIS), size = 0.7, position = position_dodge(width = .4)) +
  labs(x = 'Eccentricity (mm)', y = 'Peak speed (ms)') +
  scale_color_manual(values = c('black','grey30','grey60')) +
  facet_grid(~VIEW) + theme_classic() + 
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12)
  )


ggsave('RAD_PS_ECC.png', plot = last_plot(),  device = NULL, dpi = 300, 
       width = 7.5, height = 5, path = anaPath)

## ANOVA ##
PSECC_ANOVA <- ezANOVA(
  data = PS_ECC
  , dv = .(PS)
  , wid = .(PPT)
  , within = .(VIEW, ECC)
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
PSttest <- pairwise.t.test(PS_ECC$PS, PS_ECC$DIAGNOSIS, p.adj = 'bonf')
print(PSttest)

##### TIME TO PS #####
TPS_medians <- aggregate(TPS ~ PPT * VIEW * SIDE * DOM * POSITION * ECC * SITE * DIAGNOSIS, 
                         median, data = res)
TPS_ECC <- aggregate(TPS ~ PPT * VIEW * ECC * SITE * DIAGNOSIS, 
                         median, data = TPS_medians)
TPS_means <- aggregate(TPS ~ PPT * VIEW * SITE * DIAGNOSIS,
                        mean, data = TPS_ECC)

## plotting :)
# average across sides
TPS_means$DIAGNOSIS <- factor(TPS_means$DIAGNOSIS, levels = c('HC','MCI','AD'))

ggplot(TPS_means, aes(x = VIEW, y = TPS, colour = SITE, group = PPT)) +
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

ggsave('RAD_TPSmean.png', plot = last_plot(),  device = NULL, dpi = 300, 
       width = 6, height = 4, path = anaPath)

# by eccentricity
# summary data
TPSecc <- summarySE(TPS_medians, measurevar = 'TPS', 
                     groupvar = c('DIAGNOSIS','ECC','VIEW'), na.rm = TRUE)
TPSecc$DIAGNOSIS <- factor(TPSecc$DIAGNOSIS, levels = c('HC','MCI','AD'))

ggplot(TPSecc, aes(x = ECC, y = TPS, group = DIAGNOSIS, colour = DIAGNOSIS,
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


ggsave('RAD_TPS_ECC.png', plot = last_plot(),  device = NULL, dpi = 300, 
       width = 7.5, height = 5, path = anaPath)

## ANOVA ##
# full anova on TPS by ecc
TPSECC_ANOVA <- ezANOVA(
  data = TPS_ECC
  , dv = .(TPS)
  , wid = .(PPT)
  , within = .(VIEW, ECC)
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
TPSttest <- pairwise.t.test(TPS_ECC$TPS, TPS_ECC$DIAGNOSIS, p.adj = 'bonf')
print(TPSttest)

##### TIME AFTER PV #####
# calculating
res$TAPS <- res$MT - res$TPS

TAPS_medians <- aggregate(TAPS ~ PPT * VIEW * SIDE * DOM * POSITION * ECC * SITE * DIAGNOSIS, 
                          median, data = res)
# collapsing across side first
TAPS_ECC <- aggregate(TAPS ~ PPT * VIEW * ECC * SITE * DIAGNOSIS, 
                          median, data = TAPS_medians)
TAPS_means <- aggregate(TAPS ~ PPT * VIEW * SITE * DIAGNOSIS,
                        mean, data = TAPS_ECC)

## plotting :)
# means all PP
TAPS_means$DIAGNOSIS <- factor(TAPS_means$DIAGNOSIS, levels = c('HC','MCI','AD'))

ggplot(TAPS_means, aes(x = VIEW, y = TAPS, colour = SITE, group = PPT)) +
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

ggsave('RAD_TAPSmean.png', plot = last_plot(),  device = NULL, dpi = 300, 
       width = 6, height = 4, path = anaPath)

# by eccentricity
# summary data
TAPSecc <- summarySE(TAPS_medians, measurevar = 'TAPS', 
                     groupvar = c('DIAGNOSIS','ECC','VIEW'), na.rm = TRUE)
TAPSecc$DIAGNOSIS <- factor(TAPSecc$DIAGNOSIS, levels = c('HC','MCI','AD'))

ggplot(TAPSecc, aes(x = ECC, y = TAPS, group = DIAGNOSIS, colour = DIAGNOSIS,
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


ggsave('RAD_TAPS_ECC.png', plot = last_plot(),  device = NULL, dpi = 300, 
       width = 7.5, height = 5, path = anaPath)

## ANOVA ##
TAPSECC_ANOVA <- ezANOVA(
  data = TAPS_ECC
  , dv = .(TAPS)
  , wid = .(PPT)
  , within = .(VIEW, ECC)
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
TAPSttest <- pairwise.t.test(TAPS_ECC$TAPS, TAPS_ECC$DIAGNOSIS, p.adj = 'bonf')
print(TAPSttest)

##### NORM TIME AFTER PV #####
# calculating
res$NTAPS <- (res$MT - res$TPS)/res$MT

NTAPS_medians <- aggregate(NTAPS ~ PPT * VIEW * SIDE * DOM * POSITION * ECC * SITE * DIAGNOSIS, 
                        median, data = res)
NTAPS_ECC <- aggregate(NTAPS ~ PPT * VIEW * ECC * SITE * DIAGNOSIS, 
                       median, data = NTAPS_medians)
NTAPS_means <- aggregate(NTAPS ~ PPT * VIEW * SITE * DIAGNOSIS,
                      mean, data = NTAPS_ECC)

## plotting :)
# average across sides
NTAPS_means$DIAGNOSIS <- factor(NTAPS_means$DIAGNOSIS, levels = c('HC','MCI','AD'))

ggplot(NTAPS_means, aes(x = VIEW, y = NTAPS, colour = SITE, group = PPT)) +
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

ggsave('RAD_NTAPSmean.png', plot = last_plot(),  device = NULL, dpi = 300, 
       width = 6, height = 4, path = anaPath)

# by eccentricity
# summary data
NTAPSecc <- summarySE(NTAPS_medians, measurevar = 'NTAPS', 
                   groupvar = c('DIAGNOSIS','ECC','VIEW'), na.rm = TRUE)
NTAPSecc$DIAGNOSIS <- factor(NTAPSecc$DIAGNOSIS, levels = c('HC','MCI','AD'))

ggplot(NTAPSecc, aes(x = ECC, y = NTAPS, group = DIAGNOSIS, colour = DIAGNOSIS,
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


ggsave('RAD_NTAPS_ECC.png', plot = last_plot(),  device = NULL, dpi = 300, 
       width = 7.5, height = 5, path = anaPath)

## ANOVA ##
NTAPSECC_ANOVA <- ezANOVA(
  data = NTAPS_ECC
  , dv = .(NTAPS)
  , wid = .(PPT)
  , within = .(VIEW, ECC)
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
NTAPSttest <- pairwise.t.test(NTAPS_ECC$NTAPS, NTAPS_ECC$DIAGNOSIS, p.adj = 'bonf')
print(NTAPSttest)

##### CORRELATE ACE #####
setwd(dataPath)
patient_demos <- read.csv('patient_demographics.csv')
setwd(anaPath)
PMIdata <- read.csv('radialPMI.csv')
#extracting ACE data into seperate data-frame
ACEscores <- patient_demos[ ,c(1, 10:15)]
names(ACEscores)[1] <- 'PPT'
# merging ACE with PMI
PMIACE <- merge(PMIdata, ACEscores, by = 'PPT')
PMIACE$ACEall <- as.numeric(as.character(PMIACE$ACEall))

## aaaand correlate :)
# average PMI
PMIACE_all <- aggregate(PMI ~ PPT+DIAGNOSIS+ACEall, mean, data = PMIACE)
ggscatter(PMIACE_all, x = 'ACEall', y = 'PMI', 
          add = 'reg.line', conf.int = FALSE, add.params = list(color = "black"),
          cor.coef = TRUE, size = 1.5, cor.coef.size = 3, cor.method = 'spearman') +
  ylab('Radial PMI (deg)') + xlab('ACE score (%)') +
  theme(text = element_text(size = 10))

ggsave('RAD_PMI-ACE.png', plot = last_plot(),  device = NULL, dpi = 300, 
       width = 5, height = 5, path = anaPath)
