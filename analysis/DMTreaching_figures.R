## FIGURES FOR PUBLICATION ## 
library(readr)
library(ggplot2)
library(Rmisc)
library(gridExtra)
library(ggpubr)

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
plot_ecc$DOM <- factor(plot_ecc$DOM, labels = c('Dominant', 'Non-dominant'))
plot_ecc$POSITION <- factor(plot_ecc$POSITION)
plot_ecc$DIAGNOSIS <- factor(plot_ecc$DIAGNOSIS, levels = c('HC','MCI','AD'))
# re-labelling and ordering

# seperating, create different plots for dominant and non-dominant sides
plot_eccD <- plot_ecc[plot_ecc$DOM == 'Dominant' ,]

ggplot(plot_eccD, aes(x = POSITION, y = AEmed, colour = DIAGNOSIS, group= DIAGNOSIS)) +
  geom_point(size = 3, position = position_dodge(width = .4)) + 
  geom_errorbar(aes(ymin=AEmed-ci, ymax=AEmed+ci), 
                width=.4, position = position_dodge(width = .4)) +
  geom_line(aes(group = DIAGNOSIS), size = .8, position = position_dodge(width = .4)) +
  labs(title = 'Dominant', x = 'Eccentricity (°)', y = '', 
      element_text(size = 10)) + ylim(0,30) + 
  facet_wrap(~VIEW)  + 
  theme_classic() + theme(legend.position = 'none') -> eccD
eccD <- ggpar(eccD, palette = 'jco')

plot_eccND <- plot_ecc[plot_ecc$DOM == 'Non-dominant' ,]

ggplot(plot_eccND, aes(x = POSITION, y = AEmed, colour = DIAGNOSIS, group= DIAGNOSIS)) +
  geom_point(size = 3, position = position_dodge(width = .4)) + 
  geom_errorbar(aes(ymin=AEmed-ci, ymax=AEmed+ci), 
                width=.4, position = position_dodge(width = .4)) +
  geom_line(aes(group = DIAGNOSIS), size = .8, position = position_dodge(width = .4)) +
  labs(title = 'Non-dominant', x = 'Eccentricity (°)', y = 'Reaching error (mm)', 
       element_text(size = 10)) + ylim(0,30) + 
  facet_wrap(~VIEW)  + 
  theme_classic() + theme(legend.position = 'none') -> eccND
eccND <- ggpar(eccND, palette = 'jco')

# combining and saving
ecc <- ggarrange(eccND, eccD,
                 ncol=2, nrow=1)
ecc
ggsave('eccentricity-fig.png', plot = last_plot(), device = NULL, dpi = 300, 
       width = 8, height = 6, path = latPath)

## PLOT 2: PMI ##
# make plot data-frame
plot_PMI <- PMIdata
plot_PMI$DOM <- factor(plot_PMI$DOM, labels = c('Dominant', 'Non-dominant'))
plot_PMI$DOM <- factor(plot_PMI$DOM, levels = c('Non-dominant', 'Dominant'))
plot_PMI$DIAGNOSIS <- factor(plot_PMI$DIAGNOSIS, levels = c('HC','MCI','AD'))
plot_PMI$PPT <- factor(plot_PMI$PPT)
jitter <- position_jitter(width = 0.1, height = 0.1)

ggplot(plot_PMI, aes(x = DOM, y = PMI, colour = DIAGNOSIS, group = PPT)) + 
  geom_point(size = 4, position = position_dodge(.2)) +
  geom_line(aes(group = PPT), alpha = .5, size = .9, position = position_dodge(.2)) +
  stat_summary(aes(y = PMI, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1, size = 4, group = 1) +
  facet_wrap(~DIAGNOSIS) +
  labs(title = 'Peripheral misreaching index', x = '', y = 'Reaching error (mm)', 
       element_text(size = 14)) +
  theme_classic() + theme(legend.position = 'none', axis.text = element_text(size = 12),
                          strip.text = element_text(size = 14)) -> pPMI
pPMI <- ggpar(pPMI, palette = 'jco')
# saving
ggsave('PMI_side-fig.png', plot = last_plot(), device = NULL, dpi = 300, 
       width = 8, height = 5, path = latPath)

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
plot_ecc$DOM <- factor(plot_ecc$DOM, labels = c('Dominant', 'Non-dominant'))
plot_ecc$POSITION <- factor(plot_ecc$POSITION, labels = c('11','23','35'))
plot_ecc$DIAGNOSIS <- factor(plot_ecc$DIAGNOSIS, levels = c('HC','MCI','AD'))
# re-labelling and ordering

# seperating, create different plots for dominant and non-dominant sides
plot_eccD <- plot_ecc[plot_ecc$DOM == 'Dominant' ,]

ggplot(plot_eccD, aes(x = POSITION, y = AE, colour = DIAGNOSIS, group= DIAGNOSIS)) +
  geom_point(size = 3, position = position_dodge(width = .4)) + 
  geom_errorbar(aes(ymin=AE-ci, ymax=AE+ci), 
                width=.4, position = position_dodge(width = .4)) +
  geom_line(aes(group = DIAGNOSIS), size = .8, position = position_dodge(width = .4)) +
  labs(title = 'Dominant', x = 'Eccentricity (°)', y = '', 
       element_text(size = 10)) + ylim(-2,35) + 
  facet_wrap(~VIEW)  + 
  theme_classic() + theme(legend.position = 'none') -> eccD
eccD <- ggpar(eccD, palette = 'jco')

plot_eccND <- plot_ecc[plot_ecc$DOM == 'Non-dominant' ,]

ggplot(plot_eccND, aes(x = POSITION, y = AE, colour = DIAGNOSIS, group= DIAGNOSIS)) +
  geom_point(size = 3, position = position_dodge(width = .4)) + 
  geom_errorbar(aes(ymin=AE-ci, ymax=AE+ci), 
                width=.4, position = position_dodge(width = .4)) +
  geom_line(aes(group = DIAGNOSIS), size = .8, position = position_dodge(width = .4)) +
  labs(title = 'Non-dominant', x = 'Eccentricity (°)', y = 'Reaching error (mm)', 
       element_text(size = 10)) + 
  facet_wrap(~VIEW)  + ylim(-2,35) + 
  theme_classic() + theme(legend.position = 'none') -> eccND
eccND <- ggpar(eccND, palette = 'jco')

# combining and saving
ecc <- ggarrange(eccND, eccD,
                 ncol=2, nrow=1)
ecc
ggsave('eccentricity-fig.png', plot = last_plot(), device = NULL, dpi = 300, 
       width = 8, height = 6, path = radPath)

## PLOT 2: PMI ##
# make plot data-frame
plot_PMI <- PMIdata
plot_PMI$DOM <- factor(plot_PMI$DOM, labels = c('Dominant', 'Non-dominant'))
plot_PMI$DOM <- factor(plot_PMI$DOM, levels = c('Non-dominant', 'Dominant'))
plot_PMI$DIAGNOSIS <- factor(plot_PMI$DIAGNOSIS, levels = c('HC','MCI','AD'))
plot_PMI$PPT <- factor(plot_PMI$PPT)

ggplot(plot_PMI, aes(x = DOM, y = PMI, colour = DIAGNOSIS, group = PPT)) + 
  geom_point(size = 4, position = position_dodge(.2)) +
  geom_line(aes(group = PPT), alpha = .5, size = .9, position = position_dodge(.2)) +
  stat_summary(aes(y = PMI, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1, size = 4, group = 1) +
  facet_wrap(~DIAGNOSIS) +
  labs(title = 'Peripheral misreaching index', x = '', y = 'Reaching error (mm)', 
       element_text(size = 14)) +
  theme_classic() + theme(legend.position = 'none', axis.text = element_text(size = 12),
                          strip.text = element_text(size = 14)) -> pPMI
pPMI <- ggpar(pPMI, palette = 'jco')
# saving
ggsave('PMI_side-fig.png', plot = last_plot(), device = NULL, dpi = 300, 
       width = 8, height = 5, path = radPath)

