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
res_meds <- dcast(res_meds, PPT+GRP+SITE+DOM+POSITION+DIAGNOSIS+AGE+ED ~ VIEW)
res_meds$PMI <- res_meds$Peripheral - res_meds$Free
# make plot data-frame
plot_ecc <- summarySE(res_meds, measurevar = 'PMI', 
                      groupvar = c('DIAGNOSIS', 'DOM', 'POSITION'), na.rm = TRUE)
plot_ecc$DOM <- factor(plot_ecc$DOM, labels = c('Dominant', 'Non-dominant'))
plot_ecc$DOM <- factor(plot_ecc$DOM, levels = c('Non-dominant', 'Dominant'))
plot_ecc$POSITION <- factor(plot_ecc$POSITION)
plot_ecc$DIAGNOSIS <- factor(plot_ecc$DIAGNOSIS, levels = c('HC','MCI','AD'))
# re-labelling and ordering

# seperating, create different plots for dominant and non-dominant sides
ggplot(plot_ecc, aes(x = POSITION, y = PMI, colour = DIAGNOSIS, group= DIAGNOSIS)) +
  geom_point(size = 5, position = position_dodge(width = .5)) + 
  geom_errorbar(aes(ymin=PMI-ci, ymax=PMI+ci), 
                width=.4, position = position_dodge(width = .5)) +
  geom_line(aes(group = DIAGNOSIS), size = 1, position = position_dodge(width = .5)) +
  scale_color_manual(values = c('grey70','grey70','grey70')) +
  labs(x = 'Eccentricity (°)', y = 'PMI (mm)') + ylim(0,25) + 
  facet_wrap(~DOM) +
  theme_classic() + theme(legend.position = 'none', 
                          axis.text = element_text(size = 10),
                          axis.title = element_text(size = 12),
                          strip.text = element_text(size = 12)
  ) -> ecc

ecc

ggsave('LATeccentricity-fig.png', plot = last_plot(), device = NULL, dpi = 300, 
       width = 6, height = 3, path = latPath)

## PLOT 2: PMI ##
# make control data-frame
control_PMI <- subset(PMIdata, PMIdata$DIAGNOSIS == 'HC')
control_PMI$TSTAT <- 0
control_PMI$PVALUE <- 1
control_PMI$DEFICIT <- 0
control_PMI$BL <- 0
# get deficit data for patients
plot_PMI <- merge(PMIdata, res_cc, by = c('PPT', 'DOM', 'DIAGNOSIS', 'PMI'))
# make plot data frame
plot_PMI <- rbind(control_PMI, plot_PMI)

plot_PMI$DOM <- factor(plot_PMI$DOM, labels = c('Dom', 'Non-dom'))
plot_PMI$DOM <- factor(plot_PMI$DOM, levels = c('Non-dom', 'Dom'))
plot_PMI$DIAGNOSIS <- factor(plot_PMI$DIAGNOSIS, levels = c('HC','MCI','AD'))
plot_PMI$PPT <- factor(plot_PMI$PPT)
plot_PMI$DEFICIT <- factor(plot_PMI$DEFICIT)

ggplot(plot_PMI, aes(x = DOM, y = PMI, colour = DIAGNOSIS, group = PPT, shape = DEFICIT)) + 
  geom_point(size = 6, position = position_dodge(.2)) +
  geom_line(aes(group = PPT), alpha = .5, size = 1, position = position_dodge(.2)) +
  stat_summary(aes(y = PMI, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1.5, size = 5, group = 1) +
  scale_color_manual(values = c('grey50','goldenrod2','dodgerblue3')) +
  scale_shape_manual(values = c(16, 1)) +
  facet_wrap(~DIAGNOSIS) +
  labs(x = '', y = 'PMI (mm)') +
  theme_classic() + theme(legend.position = 'none', 
                            axis.title = element_text(size = 22),
                            axis.text = element_text(size = 20),
                            strip.text = element_text(size = 20), 
    ) -> pPMI
pPMI
# saving
ggsave('LATPMI_side-fig.png', plot = last_plot(), device = NULL, dpi = 300, 
       width = 8, height = 5, path = latPath)


## PLOT 3: scatter plot showing significant cases ##
plot_CC <- subset(PMIdata, PMIdata$GRP == 'Patient') #data-frame with just patient data
control <- subset(PMIdata, PMIdata$GRP == 'Control')
# calculate control means
control_mean <- summarySEwithin(control, measurevar = 'PMI', withinvars = c('DOM'))
# neating data-frame for plotting
plot_CC$DOM <- factor(plot_CC$DOM, labels = c('Dominant', 'Non-dominant'))
plot_CC$DOM <- factor(plot_CC$DOM, levels = c('Non-dominant', 'Dominant'))
plot_CC$DIAGNOSIS <- factor(plot_CC$DIAGNOSIS, levels = c('MCI','AD'))

# need 2 plots for this 
# DOMINANT
plot_CCD <- subset(plot_CC, plot_CC$DOM == 'Dominant')
plot_CCD$PPTindex <- factor(1:length(plot_CCD$PPT))
# getting y-intercept values for h-line
CCDy <- control_mean$PMI[1]
CCDy1 <- control_mean$sd[1]

ggplot(plot_CCD) +
  geom_hline(yintercept = CCDy, size = 1) +
  geom_hline(yintercept = CCDy-CCDy1, size = 1, linetype = 'dashed') +
  geom_hline(yintercept = CCDy+CCDy1, size = 1, linetype = 'dashed') +
  geom_point(aes(x = PPTindex, y = PMI, colour = DIAGNOSIS), size = 5) +
  scale_colour_manual(values = c('goldenrod2', 'dodgerblue3')) +
  ylim(0,50) +
  labs(title = 'Dominant', x = '', y = 'PMI (mm)') +
  theme_classic() + theme(legend.position = 'none', 
                          title = element_text(size = 22),
                          axis.title = element_text(size = 22),
                          axis.text = element_text(size = 18),
                          axis.text.x = element_text(size = 16, angle = 90, hjust = 0.5),
                          strip.text = element_text(size = 22), 
                          legend.text = element_text(size = 20),
                          axis.ticks.x = element_blank()
                          ) -> CCD

# NON-DOMINANT
plot_CCND <- subset(plot_CC, plot_CC$DOM == 'Non-dominant')
plot_CCND$PPTindex <- factor(1:length(plot_CCND$PPT))
# getting y-intercept values for h-line
CCNDy <- control_mean$PMI[2]
CCNDy1 <- control_mean$sd[2]

ggplot(plot_CCND) +
  geom_hline(yintercept = CCNDy, size = 1) +
  geom_hline(yintercept = CCNDy-CCNDy1, size = 1, linetype = 'dashed') +
  geom_hline(yintercept = CCNDy+CCNDy1, size = 1, linetype = 'dashed') +
  geom_point(aes(x = PPTindex, y = PMI, colour = DIAGNOSIS), size = 5) +
  scale_colour_manual(values = c('goldenrod2', 'dodgerblue3')) +
  labs(title = 'Non-dominant', x = '', y = 'PMI (mm)') +
  ylim(0,50) +
  theme_classic() + theme(legend.position = 'none', 
                          title = element_text(size = 22),
                          axis.title = element_text(size = 22),
                          axis.text = element_text(size = 18),
                          strip.text = element_text(size = 22), 
                          legend.text = element_text(size = 20),
                          axis.ticks.x = element_blank(),
                          axis.text.x = element_blank()
  ) -> CCND

# compiling and saving
CC <- ggarrange(CCND, CCD,
                 ncol=1, nrow=2)
CC
ggsave('LATcase-control.png', plot = last_plot(), device = NULL, dpi = 300, 
       width = 7, height = 10, path = latPath)

###### RADIAL REACHING FIGURES ######
setwd(radPath)
# loading files
res_meds <- read.csv('radial-medians-all_filtered.csv')
res_means <- read.csv('radial-means_filtered.csv')
PMIdata <- read.csv('radialPMI-filtered.csv')
res_cc <- read.csv('radial-reaching_case-control.csv')

## PLOT 1: median eccentricity ##
# make plot data-frame
# removing outer target position
res_meds$POSITION <- factor(abs(res_meds$POSITION))
res_meds <- res_meds[res_meds$POSITION != '400' ,]

res_meds <- dcast(res_meds, PPT+GRP+SITE+DOM+POSITION+DIAGNOSIS+AGE+ED ~ VIEW)
res_meds$PMI <- res_meds$Peripheral - res_meds$Free
# make plot data-frame
plot_ecc <- summarySE(res_meds, measurevar = 'PMI', 
                      groupvar = c('DIAGNOSIS', 'DOM', 'POSITION'), na.rm = TRUE)
plot_ecc$DOM <- factor(plot_ecc$DOM, labels = c('Dominant', 'Non-dominant'))
plot_ecc$DOM <- factor(plot_ecc$DOM, levels = c('Non-dominant', 'Dominant'))
plot_ecc$POSITION <- factor(plot_ecc$POSITION, labels = c('11','22','33'))
plot_ecc$DIAGNOSIS <- factor(plot_ecc$DIAGNOSIS, levels = c('HC','MCI','AD'))
# re-labelling and ordering

# seperating, create different plots for dominant and non-dominant sides
ggplot(plot_ecc, aes(x = POSITION, y = PMI, colour = DIAGNOSIS, group= DIAGNOSIS)) +
  geom_point(size = 5, position = position_dodge(width = .5)) + 
  geom_errorbar(aes(ymin=PMI-ci, ymax=PMI+ci), 
                width=.4, position = position_dodge(width = .5)) +
  geom_line(aes(group = DIAGNOSIS), size = 1, position = position_dodge(width = .5)) +
  scale_color_manual(values = c('grey50','goldenrod2','dodgerblue3')) +
  labs(x = 'Eccentricity (°)', y = 'PMI (mm)') + 
  facet_wrap(~DOM) +
  theme_classic() + theme(legend.position = 'none', 
                          axis.text = element_text(size = 18),
                          axis.title = element_text(size = 20),
                          strip.text = element_text(size = 20)
  ) -> ecc

ecc

ggsave('RADeccentricity-fig.png', plot = last_plot(), device = NULL, dpi = 300, 
       width = 8, height = 5, path = radPath)

## PLOT 2: PMI ##
# make control data-frame
control_PMI <- subset(PMIdata, PMIdata$DIAGNOSIS == 'HC')
control_PMI$TSTAT <- 0
control_PMI$PVALUE <- 1
control_PMI$DEFICIT <- 0
control_PMI$BL <- 0
# get deficit data for patients
plot_PMI <- merge(PMIdata, res_cc, by = c('PPT', 'DOM', 'DIAGNOSIS', 'PMI'))
# make plot data frame
plot_PMI <- rbind(control_PMI, plot_PMI)

plot_PMI$DOM <- factor(plot_PMI$DOM, labels = c('Dom', 'Non-dom'))
plot_PMI$DOM <- factor(plot_PMI$DOM, levels = c('Non-dom', 'Dom'))
plot_PMI$DIAGNOSIS <- factor(plot_PMI$DIAGNOSIS, levels = c('HC','MCI','AD'))
plot_PMI$PPT <- factor(plot_PMI$PPT)
plot_PMI$DEFICIT <- factor(plot_PMI$DEFICIT)

ggplot(plot_PMI, aes(x = DOM, y = PMI, colour = DIAGNOSIS, group = PPT, shape = DEFICIT)) + 
  geom_point(size = 6, position = position_dodge(.2)) +
  geom_line(aes(group = PPT), alpha = .5, size = 1, position = position_dodge(.2)) +
  stat_summary(aes(y = PMI, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1.5, size = 5, group = 1) +
  scale_color_manual(values = c('grey50','goldenrod2','dodgerblue3')) +
  scale_shape_manual(values = c(16, 1)) +
  facet_wrap(~DIAGNOSIS) +
  labs(x = '', y = 'PMI (mm)') +
  theme_classic() + theme(legend.position = 'none', 
                          axis.title = element_text(size = 22),
                          axis.text = element_text(size = 20),
                          strip.text = element_text(size = 20), 
  ) -> pPMI
pPMI
# saving
ggsave('RADPMI_side-fig.png', plot = last_plot(), device = NULL, dpi = 300, 
       width = 8, height = 5, path = radPath)

## PLOT 3: scatter plot showing significant cases ##
plot_CC <- subset(PMIdata, PMIdata$GRP == 'Patient') #data-frame with just patient data
control <- subset(PMIdata, PMIdata$GRP == 'Control')
# calculate control means
control_mean <- summarySEwithin(control, measurevar = 'PMI', withinvars = c('DOM'))
# neating data-frame for plotting
plot_CC$DOM <- factor(plot_CC$DOM, labels = c('Dominant', 'Non-dominant'))
plot_CC$DOM <- factor(plot_CC$DOM, levels = c('Non-dominant', 'Dominant'))
plot_CC$DIAGNOSIS <- factor(plot_CC$DIAGNOSIS, levels = c('MCI','AD'))

# need 2 plots for this 
# DOMINANT
plot_CCD <- subset(plot_CC, plot_CC$DOM == 'Dominant')
# adding row 12 to match non-dominant data
plot_CCD <- add_row(plot_CCD, PPT = 212)
plot_CCD <- arrange(plot_CCD, PPT)
plot_CCD$PPTindex <- factor(1:length(plot_CCD$PPT))
# getting y-intercept values for h-line
CCDy <- control_mean$PMI[1]
CCDy1 <- control_mean$sd[1]

ggplot(plot_CCD) +
  geom_hline(yintercept = CCDy, size = 1) +
  geom_hline(yintercept = CCDy-CCDy1, size = 1, linetype = 'dashed') +
  geom_hline(yintercept = CCDy+CCDy1, size = 1, linetype = 'dashed') +
  geom_point(aes(x = PPTindex, y = PMI, colour = DIAGNOSIS), size = 5) +
  scale_colour_manual(values = c('goldenrod2', 'dodgerblue3')) +
  ylim(-1,40) +
  labs(title = 'Dominant', x = '', y = 'PMI (mm)') +
  theme_classic() + theme(legend.position = 'none', 
                          title = element_text(size = 22),
                          axis.title = element_text(size = 22),
                          axis.text = element_text(size = 18),
                          axis.text.x = element_text(size = 16, angle = 90),
                          strip.text = element_text(size = 22), 
                          legend.text = element_text(size = 20),
                          axis.ticks.x = element_blank()
  ) -> CCD

# NON-DOMINANT
plot_CCND <- subset(plot_CC, plot_CC$DOM == 'Non-dominant')
plot_CCND$PPTindex <- factor(1:length(plot_CCND$PPT))
# getting y-intercept values for h-line
CCNDy <- control_mean$PMI[2]
CCNDy1 <- control_mean$sd[2]

ggplot(plot_CCND) +
  geom_hline(yintercept = CCNDy, size = 1) +
  geom_hline(yintercept = CCNDy-CCNDy1, size = 1, linetype = 'dashed') +
  geom_hline(yintercept = CCNDy+CCNDy1, size = 1, linetype = 'dashed') +
  geom_point(aes(x = PPTindex, y = PMI, colour = DIAGNOSIS), size = 5) +
  scale_colour_manual(values = c('goldenrod2', 'dodgerblue3')) +
  labs(title = 'Non-dominant', x = '', y = 'PMI (mm)') +
  ylim(-1,40) +
  theme_classic() + theme(legend.position = 'none', 
                          title = element_text(size = 22),
                          axis.title = element_text(size = 22),
                          axis.text = element_text(size = 18),
                          strip.text = element_text(size = 22), 
                          legend.text = element_text(size = 20),
                          axis.ticks.x = element_blank(),
                          axis.text.x = element_blank()
  ) -> CCND

# compiling and saving
CC <- ggarrange(CCND, CCD,
                ncol=1, nrow=2)
CC
ggsave('RADcase-control.png', plot = last_plot(), device = NULL, dpi = 300, 
       width = 7, height = 10, path = radPath)
