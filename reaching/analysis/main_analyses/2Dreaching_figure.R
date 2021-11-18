##### LATERAL REACHING ANALYSIS CODE
## Code for main analysis for Mitchell et al. ....
library(readr)
library(ggplot2)
library(Rmisc)
library(gridExtra)
library(plyr)
library(tidyverse)
library(reshape2)

#set working directory to where data is -> might need to change
anaPath <- 'S:/groups/DMT/analysis/lateral_reaching'
#anaPath <- '/Users/alex/Library/Mobile Documents/com~apple~CloudDocs/Documents/DMT/analysis/lateral_reaching'
setwd(anaPath)

# load data file
res <- read.csv('lateral-reaching_compiled.csv')
# group by x and y error for each target
LANDX <- aggregate(xerr_mm~POSITION*PPT*VIEW*SIDE*SITE*GRP*DIAGNOSIS*AGE, 
                   median, data=res)
LANDY <- aggregate(yerr_mm~POSITION*PPT*VIEW*SIDE*SITE*GRP*DIAGNOSIS*AGE, 
                   median, data=res)
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

ALLPOS_free <- ALLPOS_means[ALLPOS_means$VIEW == 'free' ,]
ALLPOS_periph <- ALLPOS_means[ALLPOS_means$VIEW == 'periph' ,]

# free
ggplot(ALLPOS_free, aes(shape = DIAGNOSIS)) +
  geom_point(aes(xerr_mm, yerr_mm)) +
  geom_point(aes(0,0), colour = 'black', shape = 1, size = 15) +
  facet_wrap(~SIDE)

# peripheral
ggplot(ALLPOS_periph, aes(shape = DIAGNOSIS)) +
  geom_point(aes(xerr_mm, yerr_mm), size = 4) +
  geom_point(aes(0,0), colour = 'black', shape = 1, size = 30) +
  facet_wrap(~SIDE)

