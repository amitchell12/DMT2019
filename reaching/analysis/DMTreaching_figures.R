## FIGURES FOR PUBLICATION ## 
library(readr)
library(ggplot2)
library(Rmisc)
library(gridExtra)
library(ggpubr)
library(tidyverse)
library(reshape2)

#set working directory to where data is
# on desktop mac
#latPath <- '/Users/Alex/Documents/DMT/analysis/lateral_reaching'
#radPath <- '/Users/Alex/Documents/DMT/analysis/radial_reaching'
# on laptop
radPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/analysis/radial_reaching'
latPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/analysis/lateral_reaching'


###### LATERAL REACHING FIGURES ######
setwd(latPath)
# loading files
res_meds <- read.csv('lateral-medians_all.csv')
res_means <- read.csv('lateral-means_all.csv')
PMIdata <- read.csv('lateralPMI_all.csv')
res_cc <- read.csv('lateral_case-control.csv')

## make PMIgrand (across sides)
PMIgrand <- aggregate(PMI ~ DIAGNOSIS*PPT*SITE*AGE*GRP, mean, data = PMIdata)

##### LAT PMI PLOTS #####
# make control data-frame
control_PMI <- subset(PMIdata, PMIdata$DIAGNOSIS == 'HC')
control_PMI$TSTAT <- 0
control_PMI$PVALUE <- 1
control_PMI$DEFICIT <- 0
control_PMI$BL <- 0
# get deficit data for patients
plot_PMI <- merge(PMIdata, res_cc, by = c('PPT','DOM','DIAGNOSIS'))
# include only relevant info
plot_PMI$PMI <- plot_PMI$PMI.x
plot_PMI <- plot_PMI[, c(1:9,12:13,21:23)]

# make plot data frame
plot_PMI <- rbind(control_PMI, plot_PMI)

plot_PMI$DOM <- factor(plot_PMI$DOM, levels = c('ND', 'D'))
plot_PMI$DIAGNOSIS <- factor(plot_PMI$DIAGNOSIS, levels = c('HC','MCI','AD'))
plot_PMI$PPT <- factor(plot_PMI$PPT)
# make column where deficit and BL cases are combined
# 1 = no deficit, 2 = borderline deficit, 3 = deficit
plot_PMI$DEFICITS <- as.numeric(plot_PMI$DEFICIT) + as.numeric(plot_PMI$BL)
plot_PMI$DEFICITS <- factor(plot_PMI$DEFICITS)

ggplot(plot_PMI, aes(x = DOM, y = PMI, colour = DIAGNOSIS, group = PPT, shape = DEFICITS)) + 
  geom_line(aes(group = PPT), alpha = .7, size = 0.7, position = position_dodge(.2)) +
  geom_point(size = 2.5, position = position_dodge(.2)) +
  stat_summary(aes(y = PMI, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1, size = 4, group = 1) +
  scale_color_manual(values = c('grey60','darkorange2','dodgerblue4')) +
  scale_shape_manual(values = c(1,18)) +
  facet_wrap(~DIAGNOSIS) +
  labs(x = 'Side', y = 'PMI (mm)') +
  theme_classic() + theme(legend.position = 'none', 
                          axis.title = element_text(size = 12),
                          axis.text = element_text(size = 10),
                          strip.text = element_text(size = 10) 
  ) -> latPMI_side
latPMI_side

## averaged across sides
PMIgrand <- PMIgrand[order(PMIgrand$PPT),]
plot_PMI$DEFICITS <- as.numeric(plot_PMI$DEFICITS)
deficit <- aggregate(DEFICITS ~ PPT*DIAGNOSIS, max, data = plot_PMI)
deficit <- deficit[order(deficit$PPT),]
PMIgrand$DEFICITS <- factor(deficit$DEFICITS)
PMIgrand$DIAGNOSIS <- factor(PMIgrand$DIAGNOSIS, levels = c('HC','MCI','AD'))

ggplot(PMIgrand, aes(x = DIAGNOSIS, y = PMI, colour = DIAGNOSIS, group = PPT, shape = DEFICITS)) + 
  geom_point(size = 4, stroke = 1, position = position_dodge(.2)) +
  stat_summary(aes(y = PMI, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 2, size = 6, group = 1) +
  scale_color_manual(values = c('grey50','darkorange2','dodgerblue4')) +
  scale_shape_manual(values = c(1,18)) + ylim(0,50) +
  labs(title = 'Lateral reaching', x = 'Diagnosis', y = 'PMI (mm)') +
  theme_classic() + theme(legend.position = 'none',
                          title = element_text(size = 18),
                          axis.title = element_text(size = 16),
                          axis.text = element_text(size = 14) 
  ) -> latPMI
latPMI

##### LAT AE PLOTS #####
# make plot data-frame
av_ecc <- summarySE(res_meds, measurevar = 'AEmed', 
                    groupvar = c('DIAGNOSIS','POSITION','VIEW'), na.rm = TRUE)

av_ecc$DIAGNOSIS <- factor(av_ecc$DIAGNOSIS, levels = c('HC','MCI','AD'))
av_ecc$POSITION <- factor(av_ecc$POSITION)

# plot 
ggplot(av_ecc, aes(x = POSITION, y = AEmed, group = DIAGNOSIS, colour = DIAGNOSIS, 
                   shape = DIAGNOSIS)) +
  geom_point(size = 4.5, position = position_dodge(width = .4)) +
  geom_errorbar(aes(ymin=AEmed-ci, ymax=AEmed+ci), 
                width=.4, position = position_dodge(width = .4)) + 
  geom_line(aes(group = DIAGNOSIS), size = 1, position = position_dodge(width = .4)) +
  scale_color_manual(values = c('grey50','darkorange2','dodgerblue4')) +
  labs(title = 'Lateral reaching', x = 'Eccentricity (°)', y = 'Absolute error (mm)') +
  facet_wrap(~VIEW) + theme_classic() + ylim(0,45) +
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        title = element_text(size = 18),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 14)
  ) -> LatAE
LatAE

##### LAT MT PLOTS #####
res <- read.csv('lateral-reaching_compiled.csv')
res$ECC <- res$POSITION
MT_medians <- aggregate(MT ~ PPT * VIEW * SIDE * ECC * SITE * DIAGNOSIS * AGE, 
                        median, data = res)
# average across side
MT_ECC <- aggregate(MT ~ PPT * VIEW * SITE * ECC * DIAGNOSIS * AGE, 
                    mean, data = MT_medians)
# mean
MT_means <- aggregate(MT ~ PPT * VIEW * SITE * DIAGNOSIS,
                      mean, data = MT_ECC)


# summary data
MTsum <- summarySE(MT_ECC, measurevar = 'MT', 
                   groupvar = c('DIAGNOSIS','ECC','VIEW'), na.rm = TRUE)
MTsum$DIAGNOSIS <- factor(MTsum$DIAGNOSIS, levels = c('HC','MCI','AD'))
MTsum$ECC <- factor(MTsum$ECC)
MTsum$VIEW <- factor(MTsum$VIEW, labels = c('Free','Peripheral'))

# plot
ggplot(MTsum, aes(x = ECC, y = MT, group = DIAGNOSIS, colour = DIAGNOSIS,
                  shape = DIAGNOSIS)) +
  geom_point(size = 4.5, position = position_dodge(width = .4)) +
  geom_errorbar(aes(ymin=MT-ci, ymax=MT+ci), 
                width=.4, position = position_dodge(width = .4)) + 
  geom_line(aes(group = DIAGNOSIS), size = 1, position = position_dodge(width = .4)) +
  labs(title = 'Lateral reaching', x = 'Eccentricity (°)', y = 'Movement time (ms)') +
  scale_color_manual(values = c('grey50','darkorange2','dodgerblue4')) +
  ylim(300,900) +
  facet_wrap(~VIEW) + theme_classic() +
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        title = element_text(size = 18),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 14)
  ) -> latMT
latMT

## MT grand means
## plotting :)
# mean all PP
MT_means$DIAGNOSIS <- factor(MT_means$DIAGNOSIS, levels = c('HC','MCI','AD'))
MT_means$PPT <- factor(MT_means$PPT)
MT_means$VIEW <- factor(MT_means$VIEW, labels = c('Free','Peripheral'))

ggplot(MT_means, aes(x = DIAGNOSIS, y = MT, group = DIAGNOSIS, colour = DIAGNOSIS)) +
  geom_boxplot(size = 0.7) +
  geom_jitter(size = 2, shape=16, position=position_jitter(0.2), alpha = .5) +
  facet_wrap(~VIEW) +
  scale_color_manual(values = c('grey50','darkorange2','dodgerblue4')) +
  ylim(250,1050) +
  labs(title = 'Lateral reaching', x = '', y = 'Movement time (ms)') +
  theme_classic() + theme(legend.position = 'none', 
                          title = element_text(size = 18),
                          axis.title = element_text(size = 16),
                          axis.text = element_text(size = 12),
                          strip.text = element_text(size = 14) 
  ) -> latMT_means
latMT_means

##### LAT RT #####
RT_medians <- aggregate(RT ~ PPT * VIEW * SIDE * ECC * SITE * DIAGNOSIS * AGE, 
                        median, data = res)
# average across side
RT_ECC <- aggregate(RT ~ PPT * VIEW * SITE * ECC * DIAGNOSIS * AGE, 
                    mean, data = RT_medians)
# mean
RT_means <- aggregate(RT ~ PPT * VIEW * SITE * DIAGNOSIS,
                      mean, data = RT_ECC)

# summary data
RTsum <- summarySE(RT_ECC, measurevar = 'RT', 
                   groupvar = c('DIAGNOSIS','ECC','VIEW'), na.rm = TRUE)
RTsum$DIAGNOSIS <- factor(RTsum$DIAGNOSIS, levels = c('HC','MCI','AD'))
RTsum$ECC <- factor(RTsum$ECC)
RTsum$VIEW <- factor(RTsum$VIEW, labels = c('Free','Peripheral'))

# plot
ggplot(RTsum, aes(x = ECC, y = RT, group = DIAGNOSIS, colour = DIAGNOSIS,
                  shape = DIAGNOSIS)) +
  geom_point(size = 4.5, position = position_dodge(width = .4)) +
  geom_errorbar(aes(ymin=RT-ci, ymax=RT+ci), 
                width=.4, position = position_dodge(width = .4)) + 
  geom_line(aes(group = DIAGNOSIS), size = 1, position = position_dodge(width = .4)) +
  labs(title = 'Lateral reaching', x = 'Eccentricity (°)', y = 'Reaction time (ms)') +
  scale_color_manual(values = c('grey50','darkorange2','dodgerblue4')) +
  ylim(300,800) +
  facet_wrap(~VIEW) + theme_classic() +
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        title = element_text(size = 18),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 14)
  ) -> latRT
latRT

# mean all PP
RT_means$DIAGNOSIS <- factor(RT_means$DIAGNOSIS, levels = c('HC','MCI','AD'))
RT_means$PPT <- factor(RT_means$PPT)
RT_means$VIEW <- factor(RT_means$VIEW, labels = c('Free','Peripheral'))

ggplot(RT_means, aes(x = DIAGNOSIS, y = RT, group = DIAGNOSIS, colour = DIAGNOSIS)) +
  geom_boxplot(size = 0.7) +
  geom_jitter(size = 2, shape=16, position=position_jitter(0.2), alpha = .5) +
  facet_wrap(~VIEW) +
  scale_color_manual(values = c('grey50','darkorange2','dodgerblue4')) +
  ylim(150,850) +
  labs(title = 'Lateral reaching', x = '', y = 'Reaction time (ms)') +
  theme_classic() + theme(legend.position = 'none', 
                          title = element_text(size = 18),
                          axis.title = element_text(size = 16),
                          axis.text = element_text(size = 12),
                          strip.text = element_text(size = 14) 
  ) -> latRT_means
latRT_means

###### RADIAL REACHING FIGURES ######
setwd(radPath)
# loading files
res_meds <- read.csv('radial-medians.csv')
res_means <- read.csv('radial-means.csv')
PMIdata <- read.csv('radialPMI.csv')
res_cc <- read.csv('radial-reaching_case-control.csv')

## calculate PMIgrand
## make PMIgrand (across sides)
PMIgrand <- aggregate(PMI ~ DIAGNOSIS*PPT*SITE*AGE*GRP, mean, data = PMIdata)

##### RAD PMI PLOTS #####
# make control data-frame
control_PMI <- subset(PMIdata, PMIdata$DIAGNOSIS == 'HC')
control_PMI$TSTAT <- 0
control_PMI$PVALUE <- 1
control_PMI$DEFICIT <- 0
control_PMI$BL <- 0
# get deficit data for patients
plot_PMI <- merge(PMIdata, res_cc, by = c('PPT','DOM','DIAGNOSIS','SITE'))
# include only relevant info
plot_PMI$PMI <- plot_PMI$PMI.x
plot_PMI <- plot_PMI[, c(1:9,23,12,13,21,22)]

# make plot data frame
plot_PMI <- rbind(control_PMI, plot_PMI)

plot_PMI$DOM <- factor(plot_PMI$DOM, levels = c('ND', 'D'))
plot_PMI$DIAGNOSIS <- factor(plot_PMI$DIAGNOSIS, levels = c('HC','MCI','AD'))
plot_PMI$PPT <- factor(plot_PMI$PPT)
# make column where deficit and BL cases are combined
# 1 = no deficit, 2 = borderline deficit, 3 = deficit
plot_PMI$DEFICITS <- as.numeric(plot_PMI$DEFICIT) + as.numeric(plot_PMI$BL)
plot_PMI$DEFICITS <- factor(plot_PMI$DEFICITS)

# actual plot
ggplot(plot_PMI, aes(x = DOM, y = PMI, colour = DIAGNOSIS, group = PPT, shape = DEFICITS)) + 
  geom_line(aes(group = PPT), alpha = .7, size = 0.7, position = position_dodge(.2)) +
  geom_point(size = 2.5, position = position_dodge(.2)) +
  stat_summary(aes(y = PMI, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1, size = 4, group = 1) +
  scale_color_manual(values = c('grey45','grey45','grey45')) +
  scale_shape_manual(values = c(1,18,18)) +
  facet_wrap(~DIAGNOSIS) +
  labs(title = 'Radial reaching', x = 'Side', y = 'PMI (mm)') +
  theme_classic() + theme(legend.position = 'none', 
                          axis.title = element_text(size = 12),
                          axis.text = element_text(size = 10),
                          strip.text = element_text(size = 10) 
  ) -> radPMI_side
radPMI_side

## average PMI across sides - combine with latPMI plot
PMIgrand <- PMIgrand[order(PMIgrand$PPT),]
plot_PMI$DEFICITS <- as.numeric(plot_PMI$DEFICITS)
deficit <- aggregate(DEFICITS ~ PPT*DIAGNOSIS, max, data = plot_PMI)
deficit <- deficit[order(deficit$PPT),]
deficit <- deficit[deficit$PPT != 310 ,]
PMIgrand$DEFICITS <- factor(deficit$DEFICITS)
PMIgrand$DIAGNOSIS <- factor(PMIgrand$DIAGNOSIS, levels = c('HC','MCI','AD'))

ggplot(PMIgrand, aes(x = DIAGNOSIS, y = PMI, colour = DIAGNOSIS, group = PPT, shape = DEFICITS)) + 
  geom_point(size = 4, stroke = 1, position = position_dodge(.2)) +
  stat_summary(aes(y = PMI, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 2, size = 6, group = 1) +
  scale_color_manual(values = c('grey50','goldenrod1','darkgreen')) +
  scale_shape_manual(values = c(1,18)) + ylim(0,50) +
  labs(title = 'Radial reaching', x = 'Diagnosis', y = 'PMI (mm)') +
  theme_classic() + theme(legend.position = 'none',
                          title = element_text(size = 18),
                          axis.title = element_text(size = 16),
                          axis.text = element_text(size = 14) 
  ) -> radPMI
radPMI

##### RAD AE PLOTS #####
# make plot data-frame
res <- read.csv('all_radial-reaching_compiled.csv')
res_medsall <- aggregate(AE ~ PPT*POSITION*VIEW*SIDE*DOM*DIAGNOSIS*GRP*SITE*AGE, 
                         median, data = res)
res_medsall$ECC <- abs(as.numeric(as.character(res_medsall$POSITION)))
res_meds_grand <- aggregate(AE ~ PPT*ECC*VIEW*DIAGNOSIS*GRP*SITE*AGE, 
                            mean, data = res_medsall)

# data-frame for plottin
av_ecc <- summarySE(res_meds_grand, measurevar = 'AE', 
                    groupvar = c('DIAGNOSIS','ECC','VIEW'), na.rm = TRUE)

av_ecc$DIAGNOSIS <- factor(av_ecc$DIAGNOSIS, levels = c('HC','MCI','AD'))
av_ecc$ECC <- factor(av_ecc$ECC)

# plot 
ggplot(av_ecc, aes(x = ECC, y = AE, group = DIAGNOSIS, colour = DIAGNOSIS, 
                   shape = DIAGNOSIS)) +
  geom_point(size = 4.5, position = position_dodge(width = .4)) +
  geom_errorbar(aes(ymin=AE-ci, ymax=AE+ci), 
                width=.4, position = position_dodge(width = .4)) + 
  geom_line(aes(group = DIAGNOSIS), size = 1, position = position_dodge(width = .4)) +
  scale_color_manual(values = c('grey50','goldenrod1','darkgreen')) +
  labs(title = 'Radial reaching', x = 'Eccentricity (mm)', y = 'Absolute error (mm)') +
  facet_wrap(~VIEW) + theme_classic() + ylim(0,45) +
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        title = element_text(size = 18),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 14)
  ) -> RadAE
RadAE

##### RAD MT PLOTS #####
res$ECC <- abs(res$POSITION)
MT_medians <- aggregate(MT ~ PPT * VIEW * SIDE * ECC * SITE * DIAGNOSIS * AGE, 
                        median, data = res)
# average across side
MT_ECC <- aggregate(MT ~ PPT * VIEW * SITE * ECC * DIAGNOSIS * AGE, 
                    mean, data = MT_medians)
# mean
MT_means <- aggregate(MT ~ PPT * VIEW * SITE * DIAGNOSIS,
                      mean, data = MT_ECC)

# summary data
MTsum <- summarySE(MT_ECC, measurevar = 'MT', 
                   groupvar = c('DIAGNOSIS','ECC','VIEW'), na.rm = TRUE)
MTsum$DIAGNOSIS <- factor(MTsum$DIAGNOSIS, levels = c('HC','MCI','AD'))
MTsum$ECC <- factor(MTsum$ECC)
MTsum$VIEW <- factor(MTsum$VIEW, labels = c('Free','Peripheral'))

# plot
ggplot(MTsum, aes(x = ECC, y = MT, group = DIAGNOSIS, colour = DIAGNOSIS,
                  shape = DIAGNOSIS)) +
  geom_point(size = 4.5, position = position_dodge(width = .4)) +
  geom_errorbar(aes(ymin=MT-ci, ymax=MT+ci), 
                width=.4, position = position_dodge(width = .4)) + 
  geom_line(aes(group = DIAGNOSIS), size = 1, position = position_dodge(width = .4)) +
  labs(title = 'Radial reaching', x = 'Eccentricity (mm)', y = 'Movement time (ms)') +
  scale_color_manual(values = c('grey50','goldenrod1','darkgreen')) +
  facet_wrap(~VIEW) + theme_classic() +
  ylim(250,900) +
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        title = element_text(size = 18),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 14)
  ) -> radMT
radMT

## grand mean MT
MT_means$DIAGNOSIS <- factor(MT_means$DIAGNOSIS, levels = c('HC','MCI','AD'))
MT_means$PPT <- factor(MT_means$PPT)
MT_means$VIEW <- factor(MT_means$VIEW, labels = c('Free','Peripheral'))

ggplot(MT_means, aes(x = DIAGNOSIS, y = MT, group = DIAGNOSIS, colour = DIAGNOSIS)) +
  geom_boxplot(size = 0.7) +
  geom_jitter(shape=16, size = 2, position=position_jitter(0.2), alpha = .5) +
  facet_wrap(~VIEW) +
  scale_color_manual(values = c('grey50','goldenrod1','darkgreen')) +
  ylim(300,1050) +
  labs(title = 'Radial reaching', x = '', y = 'Movement time (ms)') +
  theme_classic() + theme(legend.position = 'none',
                          title = element_text(size = 18),
                          axis.title = element_text(size = 16),
                          axis.text = element_text(size = 12),
                          strip.text = element_text(size = 14) 
  ) -> radMT_means
radMT_means

##### RAD RT #####
RT_medians <- aggregate(RT ~ PPT * VIEW * SIDE * ECC * SITE * DIAGNOSIS * AGE, 
                        median, data = res)
# average across side
RT_ECC <- aggregate(RT ~ PPT * VIEW * SITE * ECC * DIAGNOSIS * AGE, 
                    mean, data = RT_medians)
# mean
RT_means <- aggregate(RT ~ PPT * VIEW * SITE * DIAGNOSIS,
                      mean, data = RT_ECC)

# summary data
RTsum <- summarySE(RT_ECC, measurevar = 'RT', 
                   groupvar = c('DIAGNOSIS','ECC','VIEW'), na.rm = TRUE)
RTsum$DIAGNOSIS <- factor(RTsum$DIAGNOSIS, levels = c('HC','MCI','AD'))
RTsum$ECC <- factor(RTsum$ECC)
RTsum$VIEW <- factor(RTsum$VIEW, labels = c('Free','Peripheral'))

# plot
ggplot(RTsum, aes(x = ECC, y = RT, group = DIAGNOSIS, colour = DIAGNOSIS,
                  shape = DIAGNOSIS)) +
  geom_point(size = 4.5, position = position_dodge(width = .4)) +
  geom_errorbar(aes(ymin=RT-ci, ymax=RT+ci), 
                width=.4, position = position_dodge(width = .4)) + 
  geom_line(aes(group = DIAGNOSIS), size = 1, position = position_dodge(width = .4)) +
  labs(title = 'Radial reaching', x = 'Eccentricity (mm)', y = 'Reaction time (ms)') +
  scale_color_manual(values = c('grey50','goldenrod1','darkgreen')) +
  facet_wrap(~VIEW) + theme_classic() +
  ylim(300,800) +
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        title = element_text(size = 18),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 14)
  ) -> radRT
radRT

# mean all PP
RT_means$DIAGNOSIS <- factor(RT_means$DIAGNOSIS, levels = c('HC','MCI','AD'))
RT_means$PPT <- factor(RT_means$PPT)
RT_means$VIEW <- factor(RT_means$VIEW, labels = c('Free','Peripheral'))

ggplot(RT_means, aes(x = DIAGNOSIS, y = RT, group = DIAGNOSIS, colour = DIAGNOSIS)) +
  geom_boxplot(size = 0.7) +
  geom_jitter(size = 2, shape=16, position=position_jitter(0.2), alpha = .5) +
  facet_wrap(~VIEW) +
  scale_color_manual(values = c('grey50','goldenrod1','darkgreen')) +
  ylim(150,850) +
  labs(title = 'Radial reaching', x = '', y = 'Reaction time (ms)') +
  theme_classic() + theme(legend.position = 'none', 
                          title = element_text(size = 18),
                          axis.title = element_text(size = 16),
                          axis.text = element_text(size = 12),
                          strip.text = element_text(size = 14) 
  ) -> radRT_means
radRT_means

###### COMBINE & SAVE ######
## PMI
PMIfig <- ggarrange(latPMI, radPMI,
                    ncol=2, nrow=1,
                    widths = c(1,1),
                    hjust = -1)
PMIfig

#plotPath <- '/Users/Alex/Documents/DMT/analysis/'
plotPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/analysis/'
ggsave('PMIboth-fig.png', plot = last_plot(), device = NULL, dpi = 300, 
       width = 9.5, height = 5.5, path = plotPath)

## AE ecc
AEfig <- ggarrange(LatAE, RadAE,
                    ncol=2, nrow=1,
                    widths = c(1,1),
                    hjust = -1)
AEfig

ggsave('AEboth-fig.png', plot = last_plot(), device = NULL, dpi = 300, 
       width = 10, height = 5.5, path = plotPath)

## MT
MTfig <- ggarrange(latMT, radMT,
                   ncol=2, nrow=1,
                   widths = c(1,1),
                   hjust = -1)
MTfig

ggsave('MTboth-fig.png', plot = last_plot(), device = NULL, dpi = 300, 
       width = 10, height = 5.5, path = plotPath)
# mean
MTmean_fig <- ggarrange(latMT_means, radMT_means,
                   ncol=2, nrow=1,
                   widths = c(1,1),
                   hjust = -1)
MTmean_fig

ggsave('MTboth_means-fig.png', plot = last_plot(), device = NULL, dpi = 300, 
       width = 10, height = 5, path = plotPath)

## RT
RTfig <- ggarrange(latRT, radRT,
                   ncol=2, nrow=1,
                   widths = c(1,1),
                   hjust = -1)
RTfig

ggsave('RTboth-fig.png', plot = last_plot(), device = NULL, dpi = 300, 
       width = 10, height = 5.5, path = plotPath)

# means
RTmean_fig <- ggarrange(latRT_means, radRT_means,
                        ncol=2, nrow=1,
                        widths = c(1,1),
                        hjust = -1)
RTmean_fig

ggsave('RTboth_means-fig.png', plot = last_plot(), device = NULL, dpi = 300, 
       width = 10, height = 5, path = plotPath)
