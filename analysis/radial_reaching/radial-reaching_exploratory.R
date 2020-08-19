# EXPLORATORY ANALYSIS FOR RADIAL REACHING TASK 
## AG.Mitchell 22.04.20 ##
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
## find outliers and remove 
xclude <- read.csv('radial-reaching_outliers.csv')
res <- res[!(res$PPT %in% xclude$PPT), ]
res$ECC <- abs(res$POSITION) #adding eccentricity = absolute target position
# data-frame with only 2 target locations (middle 2)
res2 <- res[res$ECC != 100 & res$ECC != 400 ,]

##### DIRECTIONAL ERROR: median, means, PMI #####
## 2 MID TARG LOCS ##
dir_medians <- aggregate(LANDx ~ PPT * VIEW * SIDE * POSITION * SITE * GRP * DIAGNOSIS * ECC, 
                         median, data = res2)
dir_means <- aggregate(LANDx ~ PPT* VIEW * SIDE * SITE * GRP * DIAGNOSIS, 
                           mean, data = dir_medians)
# casting by task
dPMI <- dcast(dir_means, PPT+DIAGNOSIS+GRP+SIDE+SITE ~ VIEW)
dPMI$PMI <- dPMI$Peripheral - dPMI$Free
write.csv(dPMI, 'radial-reaching_dirPMI.csv', row.names = FALSE)

### PLOTTING ###
# PMI
dPMI$DIAGNOSIS <- factor(dPMI$DIAGNOSIS, levels = c('HC','MCI','AD'))
ggplot(dPMI, aes(x = SIDE, y = PMI, group = PPT, colour = SITE)) + 
  geom_point(shape = 1, size = 2, stroke = .8, position = position_dodge(width = .2)) +
  geom_line(aes(group = PPT), alpha = .5, size = .5, 
            position = position_dodge(width = .2)) +
  stat_summary(aes(y = PMI, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1, size = 3.5, group = 1) +
  scale_color_manual(values = c('dodgerblue3','grey50')) +
  labs(title = 'Mid target locs', x = 'Side', y = 'x-axis PMI (mm)', 
       element_text(size = 12)) +
  facet_wrap(~DIAGNOSIS) +
  theme_classic() + theme(legend.position = 'bottom', 
                          axis.title = element_text(size = 12),
                          axis.text = element_text(size = 10),
                          strip.text = element_text(size = 10) 
  )

ggsave('RADdPMI_plot.png', plot = last_plot(), device = NULL, dpi = 300, 
       width = 5, height = 6, path = anaPath)

# Absolute error
## plotting data frame
Err_summary <- summarySE(dir_means, measurevar = 'LANDx', 
                         groupvar = c('DIAGNOSIS', 'SIDE', 'VIEW'), na.rm = TRUE)
Err_summary$DIAGNOSIS <- factor(Err_summary$DIAGNOSIS, levels = c('HC','MCI','AD'))

ggplot(Err_summary, aes(x = SIDE, y = LANDx, colour = DIAGNOSIS, shape = DIAGNOSIS,
                        group = DIAGNOSIS)) +
  geom_point(size = 4, position = position_dodge(width = .3)) +
  geom_line(aes(group = DIAGNOSIS), size = 0.7, position = position_dodge(width = .3)) +
  geom_errorbar(aes(ymin=LANDx-ci, ymax=LANDx+ci), 
                width=.3, position = position_dodge(width = .3)) +
  scale_color_manual(values = c('black','grey30','grey60')) +
  labs(title = 'Mid target locs', 
       x = 'Side', y = 'Radial reaching error (x-axis, mm)') + 
  facet_wrap(~VIEW) +
  theme_classic() + theme(legend.position = 'bottom', 
                          legend.title = element_blank(),
                          axis.text = element_text(size = 10),
                          axis.title = element_text(size = 12),
                          strip.text = element_text(size = 12))

ggsave('RADxerr_means_plot.png', plot = last_plot(), device = NULL, dpi = 300, 
       width = 5, height = 6, path = anaPath)

### ANOVA ###
# PMI
dPMIanova <- dPMI[dPMI$PPT != 212 & dPMI$PPT != 310 
                  & dPMI$PPT != 407,] #removing participants where we only have 1 data-point

# FULL ANOVA ON FILTERED DATA
DPMI_ANOVA <- ezANOVA(
  data = dPMIanova
  , dv = .(PMI)
  , wid = .(PPT)
  , within = .(SIDE)
  , between = .(DIAGNOSIS)
  , between_covariates = .(SITE)
  , type = 3,
  return_aov = TRUE,
  detailed = TRUE
)

DPMI_ANOVA$ANOVA
DPMI_ANOVA$`Mauchly's Test for Sphericity`
DPMI_ANOVA$`Sphericity Corrections`
aovDPMI <- aovEffectSize(ezObj = DPMI_ANOVA, effectSize = "pes")
aovDispTable(aovDPMI)

# Reaching error
dECCanova <- dir_medians[dir_medians$PPT != 212 & dir_medians$PPT != 310
                         & dir_medians$PPT != 407 ,]
dECCanova$ECC <- factor(dECCanova$ECC)

# FULL ANOVA ON ECCENTRICITY DATA
DECC_ANOVA <- ezANOVA(
  data = dECCanova
  , dv = .(LANDx)
  , wid = .(PPT)
  , within = .(VIEW, SIDE, ECC)
  , between = .(DIAGNOSIS)
  , between_covariates = .(SITE)
  , type = 3,
  return_aov = TRUE,
  detailed = TRUE
)

DECC_ANOVA$ANOVA
DECC_ANOVA$`Mauchly's Test for Sphericity`
DECC_ANOVA$`Sphericity Corrections`
aovDECC <- aovEffectSize(ezObj = DECC_ANOVA, effectSize = "pes")
aovDispTable(aovDECC)


## ALL TARG LOCS ##
dir_medians_all <- aggregate(LANDx ~ PPT * VIEW * SIDE * POSITION * SITE * GRP * DIAGNOSIS * ECC, 
                         median, data = res)
dir_medians_all$PPT <- factor(dir_medians_all$PPT)

# remove outlier in P101
dir_medians_all <- dir_medians_all[!(dir_medians_all$PPT == '101' & dir_medians_all$LANDx > 50) ,]

dir_means_all <- aggregate(LANDx ~ PPT* VIEW * SIDE * SITE * GRP * DIAGNOSIS, 
                       mean, data = dir_medians_all)

# casting by task
dPMIall <- dcast(dir_means_all, PPT+DIAGNOSIS+GRP+SIDE+SITE ~ VIEW)
dPMIall$PMI <- dPMIall$Peripheral - dPMIall$Free
write.csv(dPMIall, 'radial-reaching_dirPMI_all.csv', row.names = FALSE)

### PLOTTING ###
# PMI
dPMIall$DIAGNOSIS <- factor(dPMIall$DIAGNOSIS, levels = c('HC','MCI','AD'))
ggplot(dPMIall, aes(x = SIDE, y = PMI, group = PPT, colour = SITE)) + 
  geom_point(shape = 1, size = 2, stroke = .8, position = position_dodge(width = .2)) +
  geom_line(aes(group = PPT), alpha = .5, size = .5, 
            position = position_dodge(width = .2)) +
  stat_summary(aes(y = PMI, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1, size = 3.5, group = 1) +
  scale_color_manual(values = c('dodgerblue3','grey50')) +
  labs(title = 'All target locs', x = 'Side', y = 'x-axis PMI (mm)', 
       element_text(size = 12)) +
  facet_wrap(~DIAGNOSIS) +
  theme_classic() + theme(legend.position = 'bottom', 
                          axis.title = element_text(size = 12),
                          axis.text = element_text(size = 10),
                          strip.text = element_text(size = 10) 
  )

ggsave('RADdPMI_all_plot.png', plot = last_plot(), device = NULL, dpi = 300, 
       width = 5, height = 6, path = anaPath)

# Mean absolute error
## plotting data frame
Err_summary_all <- summarySE(dir_means_all, measurevar = 'LANDx', 
                         groupvar = c('DIAGNOSIS', 'SIDE', 'VIEW'), na.rm = TRUE)
Err_summary_all$DIAGNOSIS <- factor(Err_summary_all$DIAGNOSIS, levels = c('HC','MCI','AD'))

ggplot(Err_summary_all, aes(x = SIDE, y = LANDx, colour = DIAGNOSIS, shape = DIAGNOSIS,
                        group = DIAGNOSIS)) +
  geom_point(size = 4, position = position_dodge(width = .3)) +
  geom_line(aes(group = DIAGNOSIS), size = 0.7, position = position_dodge(width = .3)) +
  geom_errorbar(aes(ymin=LANDx-ci, ymax=LANDx+ci), 
                width=.3, position = position_dodge(width = .3)) +
  scale_color_manual(values = c('black','grey30','grey60')) +
  labs(title = 'All target locs', 
       x = 'Side', y = 'Radial reaching error (x-axis, mm)') + 
  facet_wrap(~VIEW) +
  theme_classic() + theme(legend.position = 'bottom', 
                          legend.title = element_blank(),
                          axis.text = element_text(size = 10),
                          axis.title = element_text(size = 12),
                          strip.text = element_text(size = 12))

ggsave('RADxerr_allmeans_plot.png', plot = last_plot(), device = NULL, dpi = 300, 
       width = 5, height = 6, path = anaPath)

# Median absolute error
av_ecc <- summarySE(dir_medians_all, measurevar = 'LANDx', 
                    groupvar = c('DIAGNOSIS','POSITION','VIEW','SIDE'), na.rm = TRUE)
av_ecc$DIAGNOSIS <- factor(av_ecc$DIAGNOSIS, levels = c('HC','MCI','AD'))
av_ecc$POSITION <- factor(av_ecc$POSITION)

# plot 
ggplot(av_ecc, aes(x = POSITION, y = LANDx, group = DIAGNOSIS, colour = DIAGNOSIS, 
                   shape = DIAGNOSIS)) +
  geom_point(size = 3, position = position_dodge(width = .4)) +
  geom_errorbar(aes(ymin=LANDx-ci, ymax=LANDx+ci), 
                width=.4, position = position_dodge(width = .4)) + 
  geom_line(aes(group = DIAGNOSIS), size = 0.7, position = position_dodge(width = .4)) +
  geom_hline(yintercept = 0) + 
  scale_color_manual(values = c('black','grey30','grey60')) +
  labs(x = 'Eccentricity (mm)', y = 'Radial reaching error (x-axis, mm)') +
  facet_grid(~VIEW) + theme_classic() +
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12)
  )

ggsave('RADxerr_ECC_plot.png', plot = last_plot(), device = NULL, dpi = 300, 
       width = 6, height = 5, path = anaPath)

### ANOVA ###
# PMI
dPMIanova_all <- dPMIall[dPMIall$PPT != 212 & dPMIall$PPT != 310 
                  & dPMIall$PPT != 407,] #removing participants where we only have 1 data-point

# FULL ANOVA ON FILTERED DATA
ALLDPMI_ANOVA <- ezANOVA(
  data = dPMIanova_all
  , dv = .(PMI)
  , wid = .(PPT)
  , within = .(SIDE)
  , between = .(DIAGNOSIS)
  , between_covariates = .(SITE)
  , type = 3,
  return_aov = TRUE,
  detailed = TRUE
)

ALLDPMI_ANOVA$ANOVA
ALLDPMI_ANOVA$`Mauchly's Test for Sphericity`
ALLDPMI_ANOVA$`Sphericity Corrections`
aovDPMI_all <- aovEffectSize(ezObj = ALLDPMI_ANOVA, effectSize = "pes")
aovDispTable(aovDPMI_all)

# Reaching error
dECCanova_all <- dir_medians_all[dir_medians_all$PPT != 101 & dir_medians_all$PPT != 212 & 
                                   dir_medians_all$PPT != 310 & dir_medians_all$PPT != 403 ,]
dECCanova_all$ECC <- factor(dECCanova_all$ECC)

# FULL ANOVA ON ECCENTRICITY DATA
ALLDECC_ANOVA <- ezANOVA(
  data = dECCanova_all
  , dv = .(LANDx)
  , wid = .(PPT)
  , within = .(VIEW, SIDE, ECC)
  , between = .(DIAGNOSIS)
  , between_covariates = .(SITE)
  , type = 3,
  return_aov = TRUE,
  detailed = TRUE
)

ALLDECC_ANOVA$ANOVA
ALLDECC_ANOVA$`Mauchly's Test for Sphericity`
ALLDECC_ANOVA$`Sphericity Corrections`
aovDECC_all <- aovEffectSize(ezObj = ALLDECC_ANOVA, effectSize = "pes")
aovDispTable(aovDECC_all)


###### MOVEMENT TIME #######
MT_medians <- aggregate(MT ~ PPT * VIEW * SIDE * DOM * POSITION * ECC * SITE * GRP * DIAGNOSIS, 
                        median, data = res2)
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
  facet_grid(~DOM) + ylim(500,800) +
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
  facet_wrap(~VIEW) + theme_classic() + ylim(500,900) +
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
MTECCanova <- MT_medians[MT_medians$PPT != 212 & 
                           MT_medians$PPT != 310 &
                           MT_medians$PPT != 315 &
                           MT_medians$PPT != 302 ,]
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
                        median, data = res2)
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
  labs(x = '', y = 'Movement time (ms)', element_text(size = 12)) +
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
                           RT_medians$PPT != 302 ,]
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
                        median, data = res2)
PS_means <- aggregate(PS ~ PPT * VIEW * SIDE * DOM * SITE * GRP * DIAGNOSIS,
                      mean, data = PS_medians)
# average across side
PSav <- aggregate(PS ~ PPT * VIEW * SITE * GRP * DIAGNOSIS,
                  mean, data = PS_medians)

###### normalised movement time after peak speed
res$NMTPS <- (res$MT - res$TPS)/res$MT
plotNMTPS <- summarySE(data=res, measurevar = "NMTPS", 
                       groupvars = c("GRP", "POSITION", "VIEW"), na.rm = TRUE)
plotNMTPS <- na.omit(plotNMTPS)

# plot
ggplot(plotNMTPS, aes(x=POSITION, y=NMTPS, colour=GRP, group=GRP)) +
  geom_point(size=5, alpha=.5, position=position_dodge(width=.3)) +
  geom_errorbar(aes(ymin=NMTPS-ci, ymax=NMTPS+ci), width=.4, position=position_dodge(width=.3)) +
  geom_line(position=position_dodge(width=.3)) +
  facet_wrap(~VIEW) +
  theme_bw()

ggsave('normMTafterPS.png', plot = last_plot(), device = NULL, dpi = 300, 
       width = 8, height = 7, path = anaPath)

### correlate reaching error with NMTPS
# creating relevant data-frame
NMTPS_mean <- aggregate(NMTPS ~  VIEW*SIDE*PPT*GRP, mean, data = res)
NMTPS_mean <- dcast(NMTPS_mean, PPT+GRP+SIDE ~ VIEW)
# rename columns
levels(NMTPS_mean$SIDE) <- c('Left', 'Right')
levels(NMTPS_mean$GRP) <- c('Control', 'Patient') 
names(NMTPS_mean)[4] <- 'NMTPS_Free'
names(NMTPS_mean)[5] <- 'NMTPS_Periph'
NMTPS_mean$PPT <- substr(NMTPS_mean$PPT, 4, 6)
NMTPS_mean$COST <- NMTPS_mean$NMTPS_FREE - NMTPS_mean$NMTPS_PERIPH

# merge with PMI
test <- merge(PMIdata, NMTPS_mean, by = c('PPT','GRP', 'SIDE'), all = TRUE)
Periphav <- aggregate(Peripheral ~ PPT * GRP, mean, data = test)
NMTPSav <- aggregate(NMTPS_Free ~ PPT * GRP, mean, data = test)
testav <- merge(Periphav, NMTPSav, by = c('PPT', 'GRP'), all = TRUE)

# plot correlation
ggscatter(test, x = 'Peripheral', y = 'NMTPS_Free', add = 'reg.line', conf.int = TRUE,
          cor.coef = TRUE, cor.method = 'pearson') + 
  facet_grid(cols = vars(DOM), rows = vars(GRP)) + 
  ylab('Normalised MTAPS (free)') + xlab('Peripheral reaching error (deg)')

# plot correlation - averaged across sides
ggscatter(testav, x = 'Peripheral', y = 'NMTPS_Free', add = 'reg.line', conf.int = TRUE,
          cor.coef = TRUE, cor.method = 'spearman') + 
  facet_wrap(~GRP) + 
  ylab('Normalised MTAPS (free)') + xlab('Peripheral reaching error (deg)')

ggsave('normMTAPS-accuracy_correlation.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

#### time to peak speed
plotTPS <- summarySE(data=res, measurevar = "TPS", 
                     groupvars = c("GRP", "POSITION", "VIEW"), na.rm = TRUE)
plotTPS <- na.omit(plotTPS)

ggplot(plotTPS, aes(x=POSITION, y=TPS, colour=GRP, group=GRP)) +
  geom_point(size=5, alpha=.5, position=position_dodge(width=.3)) +
  geom_errorbar(aes(ymin=TPS-ci, ymax=TPS+ci), width=.4, position=position_dodge(width=.3)) +
  geom_line(position=position_dodge(width=.3)) +
  facet_wrap(~VIEW) +
  theme_bw()

