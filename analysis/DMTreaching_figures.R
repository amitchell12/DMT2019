## FIGURES FOR PUBLICATION ## 
library(readr)
library(ggplot2)
library(Rmisc)
library(gridExtra)

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
  geom_point(size = 5, alpha = .5, position = position_dodge(width = .3)) + 
  geom_errorbar(aes(ymin=AEmed-ci, ymax=AEmed+ci), 
                width=.4, position = position_dodge(width = .3)) +
  geom_line(aes(group = DIAGNOSIS), alpha = .5, size = .8, position = position_dodge(width = .3)) +
  labs(title = 'Dominant', x = 'Eccentricity (°)', y = 'Reaching error (mm)', 
      element_text(size = 10)) + ylim(0,30) + 
  facet_wrap(~VIEW)  +
  theme_bw() + theme(legend.title = element_blank()) <- eccDplot

# saving


plot_eccND <- plot_ecc[plot_ecc$DOM == 'Non-dominant' ,]

ggplot(plot_eccND, aes(x = POSITION, y = AEmed, colour = DIAGNOSIS, group= DIAGNOSIS)) +
  geom_point(size = 5, alpha = .5, position = position_dodge(width = .3)) + 
  geom_errorbar(aes(ymin=AEmed-ci, ymax=AEmed+ci), 
                width=.4, position = position_dodge(width = .3)) +
  geom_line(aes(group = DIAGNOSIS), alpha = .5, size = .8, position = position_dodge(width = .3)) +
  labs(title = 'Non-dminant', x = 'Eccentricity (°)', y = 'Reaching error (mm)', 
       element_text(size = 10)) + ylim(0,30) + 
  facet_wrap(~VIEW)  +
  theme_bw() + theme(legend.title = element_blank()) <- eccNDplot

# saving

###### RADIAL REACHING FIGURES ######