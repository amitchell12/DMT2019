## FIGURES FOR PUBLICATION ## 
library(readr)
library(ggplot2)
library(Rmisc)
library(gridExtra)
library(ggpubr)
library(tidyverse)

#set working directory to where data is
# on desktop mac
latPath <- '/Users/Alex/Documents/DMT/analysis/lateral_reaching'
radPath <- '/Users/Alex/Documents/DMT/analysis/radial_reaching'
#on pc
#latPath <- 'S:/groups/DMT/analysis/lateral_reaching'
#radPath <- 'S:/groups/DMT/analysis/radial_reaching'


###### LATERAL REACHING FIGURES ######
setwd(latPath)
# loading files
res_meds <- read.csv('lateral-medians_filtered.csv')
res_means <- read.csv('lateral-means_filtered.csv')
PMIdata <- read.csv('lateralPMI-filtered.csv')

## PLOT 1: median eccentricity ##
# make plot data-frame
plot_ecc <- summarySE(res_meds, measurevar = 'AEmed', 
                      groupvar = c('DIAGNOSIS', 'DOM', 'VIEW', 'POSITION'), na.rm = TRUE)
plot_ecc$DOM <- factor(plot_ecc$DOM, labels = c('Dom', 'Non-dom'))
plot_ecc$POSITION <- factor(plot_ecc$POSITION)
plot_ecc$DIAGNOSIS <- factor(plot_ecc$DIAGNOSIS, levels = c('HC','MCI','AD'))
# re-labelling and ordering

# seperating, create different plots for dominant and non-dominant sides
plot_eccD <- plot_ecc[plot_ecc$DOM == 'Dom' ,]

ggplot(plot_eccD, aes(x = POSITION, y = AEmed, colour = DIAGNOSIS, group= DIAGNOSIS)) +
  geom_point(size = 5, position = position_dodge(width = .5)) + 
  geom_errorbar(aes(ymin=AEmed-ci, ymax=AEmed+ci), 
                width=.4, position = position_dodge(width = .5)) +
  geom_line(aes(group = DIAGNOSIS), size = 1, position = position_dodge(width = .5)) +
  scale_color_manual(values = c('grey50','goldenrod2','dodgerblue3')) +
  labs(title = 'Dominant', x = 'Eccentricity (°)', y = 'Reaching error (mm)') + ylim(0,30) + 
  facet_wrap(~VIEW)  + 
  theme_classic() + theme(legend.position = 'none', 
                          axis.text = element_text(size = 18),
                          axis.title = element_text(size = 20),
                          strip.text = element_text(size = 20),
                          title = element_text(size = 20)
  ) -> eccD

plot_eccND <- plot_ecc[plot_ecc$DOM == 'Non-dom' ,]

ggplot(plot_eccND, aes(x = POSITION, y = AEmed, colour = DIAGNOSIS, group= DIAGNOSIS)) +
  geom_point(size = 5, position = position_dodge(width = .5)) + 
  geom_errorbar(aes(ymin=AEmed-ci, ymax=AEmed+ci), 
                width=.4, position = position_dodge(width = .5)) +
  geom_line(aes(group = DIAGNOSIS), size = 1, position = position_dodge(width = .5)) +
  scale_color_manual(values = c('grey50','goldenrod2','dodgerblue3')) +
  labs(title = 'Non-dominant', x = 'Eccentricity (°)', y = 'Reaching error (mm)') + ylim(0,30) + 
  facet_wrap(~VIEW)  + 
  theme_classic() + theme(legend.position = 'none', 
                          axis.text = element_text(size = 18),
                          axis.title = element_text(size = 20),
                          strip.text = element_text(size = 20),
                          title = element_text(size = 20)
  ) -> eccND

# combining and saving
ecc <- ggarrange(eccND, eccD,
                 ncol=2, nrow=1)
ecc
ggsave('LATeccentricity-fig.png', plot = last_plot(), device = NULL, dpi = 300, 
       width = 10, height = 6, path = latPath)

## PLOT 2: PMI ##
# make plot data-frame
plot_PMI <- PMIdata
plot_PMI$DOM <- factor(plot_PMI$DOM, labels = c('Dom', 'Non-dom'))
plot_PMI$DOM <- factor(plot_PMI$DOM, levels = c('Non-dom', 'Dom'))
plot_PMI$DIAGNOSIS <- factor(plot_PMI$DIAGNOSIS, levels = c('HC','MCI','AD'))
plot_PMI$PPT <- factor(plot_PMI$PPT)
jitter <- position_jitter(width = 0.1, height = 0.1)

ggplot(plot_PMI, aes(x = DOM, y = PMI, colour = DIAGNOSIS, group = PPT)) + 
  geom_point(size = 5, position = position_dodge(.2)) +
  geom_line(aes(group = PPT), alpha = .5, size = 1, position = position_dodge(.2)) +
  stat_summary(aes(y = PMI, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1.5, size = 4.5, group = 1) +
  scale_color_manual(values = c('grey50','goldenrod2','dodgerblue3')) +
  facet_wrap(~DIAGNOSIS) +
  labs(x = '', y = 'PMI (mm)') +
  theme_classic() + theme(legend.position = 'right', legend.title = element_blank(),
                            axis.title = element_text(size = 22),
                            axis.text = element_text(size = 18),
                            strip.text = element_text(size = 22), 
                            legend.text = element_text(size = 20)) -> pPMI
# saving
ggsave('LATPMI_side-fig.png', plot = last_plot(), device = NULL, dpi = 300, 
       width = 10, height = 6, path = latPath)


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

## PLOT 1: median eccentricity ##
# make plot data-frame
# removing outer target position
res_meds$POSITION <- factor(abs(res_meds$POSITION))
res_meds <- res_meds[res_meds$POSITION != '400' ,]

plot_ecc <- summarySE(res_meds, measurevar = 'AE', 
                      groupvar = c('DIAGNOSIS', 'DOM', 'VIEW', 'POSITION'), na.rm = FALSE)
plot_ecc$DOM <- factor(plot_ecc$DOM, labels = c('Dom', 'Non-dom'))
plot_ecc$POSITION <- factor(plot_ecc$POSITION, labels = c('11','22','33'))
plot_ecc$DIAGNOSIS <- factor(plot_ecc$DIAGNOSIS, levels = c('HC','MCI','AD'))
# re-labelling and ordering

# seperating, create different plots for dominant and non-dominant sides
# seperating, create different plots for dominant and non-dominant sides
plot_eccD <- plot_ecc[plot_ecc$DOM == 'Dom' ,]

ggplot(plot_eccD, aes(x = POSITION, y = AE, colour = DIAGNOSIS, group= DIAGNOSIS)) +
  geom_point(size = 5, position = position_dodge(width = .5)) + 
  geom_errorbar(aes(ymin=AE-ci, ymax=AE+ci), 
                width=.4, position = position_dodge(width = .5)) +
  geom_line(aes(group = DIAGNOSIS), size = 1, position = position_dodge(width = .5)) +
  scale_color_manual(values = c('grey50','goldenrod2','dodgerblue3')) +
  labs(title = 'Dominant', x = 'Eccentricity (°)', y = 'Reaching error (mm)') + ylim(-10,35) + 
  facet_wrap(~VIEW)  + 
  theme_classic() + theme(legend.position = 'none', 
                          axis.text = element_text(size = 18),
                          axis.title = element_text(size = 20),
                          strip.text = element_text(size = 20),
                          title = element_text(size = 20)
  ) -> eccD

plot_eccND <- plot_ecc[plot_ecc$DOM == 'Non-dom' ,]

ggplot(plot_eccND, aes(x = POSITION, y = AE, colour = DIAGNOSIS, group= DIAGNOSIS)) +
  geom_point(size = 5, position = position_dodge(width = .5)) + 
  geom_errorbar(aes(ymin=AE-ci, ymax=AE+ci), 
                width=.4, position = position_dodge(width = .5)) +
  geom_line(aes(group = DIAGNOSIS), size = 1, position = position_dodge(width = .5)) +
  scale_color_manual(values = c('grey50','goldenrod2','dodgerblue3')) +
  labs(title = 'Non-dominant', x = 'Eccentricity (°)', y = 'Reaching error (mm)') + ylim(-10,35) + 
  facet_wrap(~VIEW)  + 
  theme_classic() + theme(legend.position = 'none', 
                          axis.text = element_text(size = 18),
                          axis.title = element_text(size = 20),
                          strip.text = element_text(size = 20),
                          title = element_text(size = 20)
  ) -> eccND

# combining and saving
ecc <- ggarrange(eccND, eccD,
                 ncol=2, nrow=1)
ecc
ggsave('RADeccentricity-fig.png', plot = last_plot(), device = NULL, dpi = 300, 
       width = 10, height = 6, path = radPath)

## PLOT 2: PMI ##
# make plot data-frame
plot_PMI <- PMIdata
plot_PMI$DOM <- factor(plot_PMI$DOM, labels = c('Dom', 'Non-dom'))
plot_PMI$DOM <- factor(plot_PMI$DOM, levels = c('Non-dom', 'Dom'))
plot_PMI$DIAGNOSIS <- factor(plot_PMI$DIAGNOSIS, levels = c('HC','MCI','AD'))
plot_PMI$PPT <- factor(plot_PMI$PPT)

ggplot(plot_PMI, aes(x = DOM, y = PMI, colour = DIAGNOSIS, group = PPT)) + 
  geom_point(size = 5, position = position_dodge(.2)) +
  geom_line(aes(group = PPT), alpha = .5, size = 1, position = position_dodge(.2)) +
  stat_summary(aes(y = PMI, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1.5, size = 4.5, group = 1) +
  scale_color_manual(values = c('grey50','goldenrod2','dodgerblue3')) +
  facet_wrap(~DIAGNOSIS) +
  labs(x = '', y = 'PMI (mm)') +
  theme_classic() + theme(legend.position = 'right', legend.title = element_blank(),
                          axis.title = element_text(size = 22),
                          axis.text = element_text(size = 18),
                          strip.text = element_text(size = 22), 
                          legend.text = element_text(size = 20)) -> pPMI
# saving
ggsave('RADPMI_side-fig.png', plot = last_plot(), device = NULL, dpi = 300, 
       width = 10, height = 6, path = radPath)

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
