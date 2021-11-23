##### LATERAL REACHING ANALYSIS CODE
## Code for main analysis for Mitchell et al. ....
library(readr)
library(ggplot2)
library(Rmisc)
library(gridExtra)
library(plyr)
library(tidyverse)
library(ggpubr)
library(reshape2)

#set working directory to where data is -> might need to change
anaPath <- 'S:/groups/DMT/analysis/lateral_reaching'
#anaPath <- '/Users/alex/Library/Mobile Documents/com~apple~CloudDocs/Documents/DMT/analysis/lateral_reaching'
setwd(anaPath)

# load data file
res <- read.csv('lateral-reaching_compiled.csv')
# group by x and y error for each target - first by eccentricity (median)
LANDX <- aggregate(xerr_mm~POSITION*PPT*VIEW*SIDE*SITE*GRP*DIAGNOSIS*AGE, 
                   median, data=res)
LANDY <- aggregate(yerr_mm~POSITION*PPT*VIEW*SIDE*SITE*GRP*DIAGNOSIS*AGE, 
                   median, data=res)
# then by side (mean)
LANDX <- aggregate(xerr_mm~PPT*VIEW*SIDE*SITE*GRP*DIAGNOSIS*AGE, 
                   mean, data=LANDX)
LANDY <- aggregate(yerr_mm~PPT*VIEW*SIDE*SITE*GRP*DIAGNOSIS*AGE, 
                   mean, data=LANDY)

ALLPOS <- merge(LANDX, LANDY)

# mean error for each target per position x & y - not absolute
Xpos_means <- summarySEwithin(data = ALLPOS, measurevar = 'xerr_mm', 
                              betweenvars = 'DIAGNOSIS', withinvars = c('VIEW','SIDE'))
names(Xpos_means)[8] <- 'CIX'
Ypos_means <- summarySEwithin(data = ALLPOS, measurevar = 'yerr_mm', 
                              betweenvars = 'DIAGNOSIS', withinvars = c('VIEW','SIDE'))
names(Ypos_means)[8] <- 'CIY'

ALLPOS_means <- merge(Xpos_means, Ypos_means, 
                      by = c('DIAGNOSIS', 'VIEW', 'SIDE', 'N'))
ALLPOS_means$DIAGNOSIS <- factor(ALLPOS_means$DIAGNOSIS, levels = c('HC','MCI','AD'))
ALLPOS_free <- ALLPOS_means[ALLPOS_means$VIEW == 'free' ,]
ALLPOS_periph <- ALLPOS_means[ALLPOS_means$VIEW == 'periph' ,]

# free -left
ggplot(ALLPOS_free[ALLPOS_free$SIDE=='left' ,], 
       aes(shape = DIAGNOSIS, colour = DIAGNOSIS)) +
  geom_point(aes(xerr_mm, yerr_mm), size = 4.5) +
  geom_errorbar(aes(x = xerr_mm, ymin=yerr_mm-se.y, ymax=yerr_mm+se.y), 
                width = .17, size = .75) +
  geom_errorbarh(aes(xmin=xerr_mm-se.x, xmax=xerr_mm+se.x, y = yerr_mm), 
                 height = .1, size = .75) +
  geom_hline(yintercept = 0, size = .75, linetype = 'dashed') + 
  geom_vline(xintercept = 0, size = .75, linetype = 'dashed') +
  geom_point(aes(0,0), fill = 'white', shape = 21, colour = 'black', 
             size = 10, stroke = 1.5)  +
  scale_colour_manual(values = c('black','grey30','grey60')) +
  ylim(-1,3) + xlim(-1,3) +
  labs(x = 'Error (x-axis, mm)', y = 'Error (y-axis, mm)', title = 'Left Free') +
  theme_classic() + theme(legend.position = 'none') -> LFP

# free -right
ggplot(ALLPOS_free[ALLPOS_free$SIDE=='right' ,], 
       aes(shape = DIAGNOSIS, colour = DIAGNOSIS)) +
  geom_point(aes(xerr_mm, yerr_mm), size = 4.5) +
  geom_errorbar(aes(x = xerr_mm, ymin=yerr_mm-se.y, ymax=yerr_mm+se.y), 
                width = .17, size = .75) +
  geom_errorbarh(aes(xmin=xerr_mm-se.x, xmax=xerr_mm+se.x, y = yerr_mm), 
                 height = .1, size = .75) +
  geom_hline(yintercept = 0, size = .75, linetype = 'dashed') + 
  geom_vline(xintercept = 0, size = .75, linetype = 'dashed') +
  geom_point(aes(0,0), fill = 'white', shape = 21, colour = 'black', 
             size = 10, stroke = 1.5)  +
  scale_colour_manual(values = c('black','grey30','grey60')) +
  ylim(-1,3) + xlim(-1,3) +
  labs(x = 'Error (x-axis, mm)', y = 'Error (y-axis, mm)', title = 'Right Free') +
  theme_classic() + theme(legend.position = c(.8,.85)) -> RFP

# peripheral -left
ggplot(ALLPOS_periph[ALLPOS_periph$SIDE=='left' ,], 
       aes(shape = DIAGNOSIS, colour = DIAGNOSIS)) +
  geom_point(aes(xerr_mm, yerr_mm), size = 4.5) +
  geom_errorbar(aes(x = xerr_mm, ymin=yerr_mm-se.y, ymax=yerr_mm+se.y), 
                width = .17, size = .75) +
  geom_errorbarh(aes(xmin=xerr_mm-se.x, xmax=xerr_mm+se.x, y = yerr_mm), 
                 height = .1, size = .75) +
  geom_hline(yintercept = 0, size = .75, linetype = 'dashed') + 
  geom_vline(xintercept = 0, size = .75, linetype = 'dashed') +
  geom_point(aes(0,0), fill = 'white', shape = 21, colour = 'black', 
             size = 10, stroke = 1.5)  +
  scale_colour_manual(values = c('black','grey30','grey60')) +
  ylim(-6,1) + xlim(-1,12) +
  labs(x = 'Error (x-axis, mm)', y = 'Error (y-axis, mm)', title = 'Left Peripheral') +
  theme_classic() + theme(legend.position = 'none') -> LPP

# peripheral -right
ggplot(ALLPOS_periph[ALLPOS_periph$SIDE=='right' ,], 
       aes(shape = DIAGNOSIS, colour = DIAGNOSIS)) +
  geom_point(aes(xerr_mm, yerr_mm), size = 4.5) +
  geom_errorbar(aes(x = xerr_mm, ymin=yerr_mm-se.y, ymax=yerr_mm+se.y), 
                width = .17, size = .75) +
  geom_errorbarh(aes(xmin=xerr_mm-se.x, xmax=xerr_mm+se.x, y = yerr_mm), 
                 height = .1, size = .75) +
  geom_hline(yintercept = 0, size = .75, linetype = 'dashed') + 
  geom_vline(xintercept = 0, size = .75, linetype = 'dashed') +
  geom_point(aes(0,0), fill = 'white', shape = 21, colour = 'black', 
             size = 10, stroke = 1.5)  +
  scale_colour_manual(values = c('black','grey30','grey60')) +
  ylim(-6,1) + xlim(-12,1) +
  labs(x = 'Error (x-axis, mm)', y = 'Error (y-axis, mm)', title = 'Right Peripheral') +
  theme_classic() + theme(legend.position = 'none') -> RPP

## combining all
TWOD <- ggarrange(LFP,RFP,LPP,RPP,
                  ncol = 2, nrow = 2,
                  widths = c(1,1),
                  hjust = .1)
TWOD
# saving
ggsave('2DREACH.png', plot = last_plot(),  device = NULL, dpi = 300, 
       width = 8, height = 8, path = anaPath)

