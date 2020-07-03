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
latPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/analysis/lateral_reaching'
radPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/analysis/radial_reaching'
#on pc
#latPath <- 'S:/groups/DMT/analysis/lateral_reaching'
#radPath <- 'S:/groups/DMT/analysis/radial_reaching'


###### LATERAL REACHING FIGURES ######
setwd(latPath)
# loading files
res_meds <- read.csv('lateral-medians_filtered.csv')
res_means <- read.csv('lateral-means_filtered.csv')
PMIdata <- read.csv('lateralPMI-filtered.csv')
res_cc <- read.csv('lateral-reaching_case-control.csv')

## PLOT 1: median eccentricity ##
# calculate PMI for each eccentricity
dom_meds <- res_meds[res_meds$DOM == 'D' ,]
ndom_meds <- res_meds[res_meds$DOM == 'ND' ,]

## Non-dominant
# make plot data-frame
plot_NDecc <- summarySE(ndom_meds, measurevar = 'AEmed', 
                        groupvar = c('DIAGNOSIS','POSITION','VIEW'), na.rm = TRUE)
plot_NDecc$POSITION <- factor(plot_NDecc$POSITION)
plot_NDecc$DIAGNOSIS <- factor(plot_NDecc$DIAGNOSIS, levels = c('HC','MCI','AD'))


# plot
ggplot(plot_NDecc, aes(x = POSITION, y = AEmed, shape = DIAGNOSIS, colour = DIAGNOSIS,
                       group= DIAGNOSIS)) +
  geom_point(size = 4, position = position_dodge(width = .3)) + 
  geom_errorbar(aes(ymin=AEmed-ci, ymax=AEmed+ci), 
                width=.4, position = position_dodge(width = .3)) +
  geom_line(aes(group = DIAGNOSIS), size = 0.7, position = position_dodge(width = .3)) +
  scale_shape_manual(values = c(1, 16, 16)) + ylim(0,30) +
  scale_color_manual(values = c('black','grey60','grey20')) +
  labs(title = 'Non-dominant', x = 'Eccentricity (°)', y = 'PMI (mm)') + 
  facet_wrap(~VIEW) +
  theme_classic() + theme(legend.position = 'none', 
                          legend.title = element_blank(),
                          axis.text = element_text(size = 10),
                          axis.title = element_text(size = 12),
                          strip.text = element_text(size = 12)
  ) -> NDecc

NDecc

## Dominant
# make plot data-frame
plot_Decc <- summarySE(dom_meds, measurevar = 'AEmed', 
                        groupvar = c('DIAGNOSIS','POSITION','VIEW'), na.rm = TRUE)
plot_Decc$POSITION <- factor(plot_Decc$POSITION)
plot_Decc$DIAGNOSIS <- factor(plot_Decc$DIAGNOSIS, levels = c('HC','MCI','AD'))


# plot
ggplot(plot_Decc, aes(x = POSITION, y = AEmed, shape = DIAGNOSIS, colour = DIAGNOSIS,
                       group= DIAGNOSIS)) +
  geom_point(size = 4, position = position_dodge(width = .3)) + 
  geom_errorbar(aes(ymin=AEmed-ci, ymax=AEmed+ci), 
                width=.4, position = position_dodge(width = .3)) +
  geom_line(aes(group = DIAGNOSIS), size = 0.7, position = position_dodge(width = .3)) +
  scale_shape_manual(values = c(1, 16, 16)) + ylim(0,30) +
  scale_color_manual(values = c('black','grey60','grey20')) +
  labs(title = 'Dominant', x = 'Eccentricity (°)', y = 'PMI (mm)') + 
  facet_wrap(~VIEW) +
  theme_classic() + theme(legend.position = 'none', 
                          legend.title = element_blank(),
                          axis.text = element_text(size = 10),
                          axis.title = element_text(size = 12),
                          strip.text = element_text(size = 12)
  ) -> Decc

Decc

ecc <- ggarrange(NDecc, Decc,
                 ncol=1, nrow=2)
ecc

ggsave('LATeccentricity-fig.png', plot = last_plot(), device = NULL, dpi = 300, 
       width = 5, height = 8, path = latPath)

## PLOT 2: PMI ##
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
plot_PMI <- plot_PMI[, c(1:10,13:14,22:24)]

# make plot data frame
plot_PMI <- rbind(control_PMI, plot_PMI)

plot_PMI$DOM <- factor(plot_PMI$DOM, labels = c('D', 'ND'))
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
  scale_color_manual(values = c('grey45','grey45','grey45')) +
  scale_shape_manual(values = c(1,18,16)) +
  facet_wrap(~DIAGNOSIS) +
  labs(x = 'Side', y = 'PMI (mm)') +
  theme_classic() + theme(legend.position = 'none', 
                            axis.title = element_text(size = 12),
                            axis.text = element_text(size = 10),
                            strip.text = element_text(size = 10) 
    ) -> pPMI
pPMI

## PLOT 3: average PMI across sides - combine with PLOT 2
PMIav_plot <- aggregate(PMI ~ PPT*DIAGNOSIS*AGE*ED*SITE, mean, data = plot_PMI)
PMIav_plot <- PMIav_plot[order(PMIav_plot$PPT),]
plot_PMI$DEFICITS <- as.numeric(plot_PMI$DEFICITS)
deficit <- aggregate(DEFICITS ~ PPT*DIAGNOSIS, max, data = plot_PMI)
deficit <- deficit[order(deficit$PPT),]
PMIav_plot$DEFICITS <- factor(deficit$DEFICITS)

ggplot(PMIav_plot, aes(x = DIAGNOSIS, y = PMI, colour = DIAGNOSIS, group = PPT, shape = DEFICITS)) + 
  geom_point(size = 2.5, position = position_dodge(.2)) +
  stat_summary(aes(y = PMI, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1, size = 4, group = 1) +
  scale_color_manual(values = c('grey45','grey45','grey45')) +
  scale_shape_manual(values = c(1,18,16)) +
  labs(x = 'Diagnosis', y = '') +
  theme_classic() + theme(legend.position = 'none', 
                          axis.title = element_text(size = 12),
                          axis.text = element_text(size = 10),
                          strip.text = element_text(size = 10) 
  ) -> avPMI
avPMI
  
PMIfig <- ggarrange(pPMI, avPMI,
                ncol=2, nrow=1,
                widths = c(1.5,1),
                labels = c('a','b'),
                hjust = -1)
PMIfig
ggsave('LATPMI-fig.png', plot = last_plot(), device = NULL, dpi = 300, 
       width = 8, height = 4, path = latPath)

###### RADIAL REACHING FIGURES ######
setwd(radPath)
# loading files
res_meds <- read.csv('radial-medians-all_filtered.csv')
res_means <- read.csv('radial-means_filtered.csv')
PMIdata <- read.csv('radialPMI-filtered.csv')
res_cc <- read.csv('radial-reaching_case-control.csv')

## PLOT 1: median eccentricity ##
# getting eccentricity info
res_meds$ECC <- factor(abs(res_meds$POSITION))
res_meds <- res_meds[res_meds$ECC != '400' ,]
# calculate PMI for each eccentricity
dom_meds <- res_meds[res_meds$DOM == 'D' ,]
ndom_meds <- res_meds[res_meds$DOM == 'ND' ,]

## Non-dominant
# make plot data-frame
plot_NDecc <- summarySE(ndom_meds, measurevar = 'AE', 
                      groupvar = c('DIAGNOSIS','ECC','VIEW'), na.rm = TRUE)
plot_NDecc$ECC <- factor(plot_NDecc$ECC,labels = c('11','22','33'))
plot_NDecc$DIAGNOSIS <- factor(plot_NDecc$DIAGNOSIS, levels = c('HC','MCI','AD'))
# re-labelling and ordering

# plot
ggplot(plot_NDecc, aes(x = ECC, y = AE, shape = DIAGNOSIS, colour = DIAGNOSIS,
                     group= DIAGNOSIS)) +
  geom_point(size = 4, position = position_dodge(width = .3)) + 
  geom_errorbar(aes(ymin=AE-ci, ymax=AE+ci), 
                width=.4, position = position_dodge(width = .3)) +
  geom_line(aes(group = DIAGNOSIS), size = 0.7, position = position_dodge(width = .3)) +
  scale_shape_manual(values = c(1, 16, 16)) + ylim(0,30) +
  scale_color_manual(values = c('black','grey60','grey20')) +
  labs(title = 'Non-dominant', x = 'Eccentricity (°)', y = 'PMI (mm)') + 
  facet_wrap(~VIEW) +
  theme_classic() + theme(legend.position = 'none', 
                          legend.title = element_blank(),
                          axis.text = element_text(size = 10),
                          axis.title = element_text(size = 12),
                          strip.text = element_text(size = 12)
  ) -> NDecc

NDecc

## Dominant
# make plot data-frame
plot_Decc <- summarySE(dom_meds, measurevar = 'AE', 
                        groupvar = c('DIAGNOSIS','ECC','VIEW'), na.rm = TRUE)
plot_Decc$ECC <- factor(plot_Decc$ECC,labels = c('11','22','33'))
plot_Decc$DIAGNOSIS <- factor(plot_Decc$DIAGNOSIS, levels = c('HC','MCI','AD'))
# re-labelling and ordering

# plot
ggplot(plot_Decc, aes(x = ECC, y = AE, shape = DIAGNOSIS, colour = DIAGNOSIS,
                       group= DIAGNOSIS)) +
  geom_point(size = 4, position = position_dodge(width = .3)) + 
  geom_errorbar(aes(ymin=AE-ci, ymax=AE+ci), 
                width=.4, position = position_dodge(width = .3)) +
  geom_line(aes(group = DIAGNOSIS), size = 0.7, position = position_dodge(width = .3)) +
  scale_shape_manual(values = c(1, 16, 16)) + ylim(0,30) +
  scale_color_manual(values = c('black','grey60','grey20')) +
  labs(title = 'Dominant', x = 'Eccentricity (°)', y = 'PMI (mm)') + 
  facet_wrap(~VIEW) +
  theme_classic() + theme(legend.position = 'none', 
                          legend.title = element_blank(),
                          axis.text = element_text(size = 10),
                          axis.title = element_text(size = 12),
                          strip.text = element_text(size = 12)
  ) -> Decc

Decc

ecc <- ggarrange(NDecc, Decc,
          ncol=1, nrow=2)
ecc

ggsave('RADeccentricity-fig.png', plot = last_plot(), device = NULL, dpi = 300, 
       width = 5, height = 8, path = radPath)

## PLOT 2: PMI ##
# make control data-frame
control_PMI <- subset(PMIdata, PMIdata$DIAGNOSIS == 'HC')
control_PMI$TSTAT <- 0
control_PMI$PVALUE <- 1
control_PMI$DEFICIT <- 0
control_PMI$BL <- 0
# get deficit data for patients
plot_PMI <- merge(PMIdata, res_cc, by = c('PPT', 'DOM', 'DIAGNOSIS', 'PMI','SITE'))
# include only relevant info
plot_PMI <- plot_PMI[, c(1:13,21:22)]

# make plot data frame
plot_PMI <- rbind(control_PMI, plot_PMI)

plot_PMI$DOM <- factor(plot_PMI$DOM, labels = c('Dom', 'Non-dom'))
plot_PMI$DOM <- factor(plot_PMI$DOM, levels = c('Non-dom', 'Dom'))
plot_PMI$DIAGNOSIS <- factor(plot_PMI$DIAGNOSIS, levels = c('HC','MCI','AD'))
plot_PMI$PPT <- factor(plot_PMI$PPT)
plot_PMI$DEFICIT <- factor(plot_PMI$DEFICIT)

ggplot(plot_PMI, aes(x = DOM, y = PMI, colour = SITE, group = PPT, shape = DEFICIT)) + 
  geom_line(aes(group = PPT), alpha = .7, size = 0.7, position = position_dodge(.2)) +
  geom_point(size = 2.5, position = position_dodge(.2)) +
  stat_summary(aes(y = PMI, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1, size = 4, group = 1) +
  scale_color_manual(values = c('black','grey45')) +
  scale_shape_manual(values = c(1,18,16)) +
  facet_wrap(~DIAGNOSIS) +
  labs(x = 'Side', y = 'PMI (mm)') +
  theme_classic() + theme(legend.position = 'none', 
                          axis.title = element_text(size = 12),
                          axis.text = element_text(size = 10),
                          strip.text = element_text(size = 10) 
  ) -> pPMI
pPMI

## PLOT 3: average PMI across sides - combine with PLOT 2
PMIav_plot <- aggregate(PMI ~ PPT*DIAGNOSIS*AGE*ED*SITE, mean, data = plot_PMI)
PMIav_plot <- PMIav_plot[order(PMIav_plot$PPT),]
plot_PMI$DEFICIT <- as.numeric(plot_PMI$DEFICIT)
deficit <- aggregate(DEFICIT ~ PPT*DIAGNOSIS, max, data = plot_PMI)
deficit <- deficit[order(deficit$PPT),]
PMIav_plot$DEFICIT <- factor(deficit$DEFICIT)

ggplot(PMIav_plot, aes(x = DIAGNOSIS, y = PMI, colour = SITE, group = PPT, shape = DEFICIT)) + 
  geom_point(size = 2.5, position = position_dodge(.2)) +
  stat_summary(aes(y = PMI, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1, size = 4, group = 1) +
  scale_color_manual(values = c('black','grey45')) +
  scale_shape_manual(values = c(1,18,16)) +
  labs(x = 'Diagnosis', y = '') +
  theme_classic() + theme(legend.position = 'none', 
                          axis.title = element_text(size = 12),
                          axis.text = element_text(size = 10),
                          strip.text = element_text(size = 10) 
  ) -> avPMI
avPMI

PMIfig <- ggarrange(pPMI, avPMI,
                    ncol=2, nrow=1,
                    widths = c(1.5,1),
                    labels = c('a','b'),
                    hjust = -1)
PMIfig

ggsave('RADPMI-fig.png', plot = last_plot(), device = NULL, dpi = 300, 
       width = 8, height = 4, path = radPath)
