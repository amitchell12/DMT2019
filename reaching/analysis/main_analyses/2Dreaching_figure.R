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
#anaPath <- 'S:/groups/DMT/analysis/lateral_reaching'
anaPath <- '/Users/alex/Library/Mobile Documents/com~apple~CloudDocs/Documents/DMT/analysis/lateral_reaching'
setwd(anaPath)

###### SHAPING DATA #####
# load data file
res <- read.csv('lateral-reaching_compiled.csv')

# screen information
x = 310
y = 175
pixels_perdeg = 43.76
deg_perpix = 1/pixels_perdeg
x_res = 1920
y_res = 1080
sd = 40

pixPer_mm_x = x_res/x;
pixPer_mm_y = y_res/y;
pixPer_mm = (pixPer_mm_x+pixPer_mm_y)/2;
mm_perPix = 1/pixPer_mm;

# convert LAND/TARG from pix to mm
res$LANDx_mm <- res$LANDx*mm_perPix
res$LANDy_mm <- res$LANDy*mm_perPix
res$TARGx_mm <- res$TARGx*mm_perPix
res$TARGy_mm <- res$TARGy*mm_perPix

# grouping
res$ECCx <- factor(cut(abs(res$TARGx_mm), 3), labels = c('49', '83', '118'))
#res$ECCy <- factor(cut(abs(res$TARGy_mm), 3), labels = c('-30', '0', '30'))

TARGx <- summarySEwithin(data = res, measurevar = 'TARGx_mm', withinvars = c('ECCx','SIDE'))

# group by x and y error for each target - first by eccentricity (median)
LANDX <- aggregate(LANDx_mm~POSITION*ECCx*PPT*VIEW*SIDE*SITE*GRP*DIAGNOSIS*AGE, 
                   median, data=res)
LANDY <- aggregate(LANDy_mm~POSITION*ECCx*PPT*VIEW*SIDE*SITE*GRP*DIAGNOSIS*AGE, 
                   median, data=res)

ALLPOS <- merge(LANDX, LANDY)

# isolate peripheral condition
ALLPERIPH <- ALLPOS[ALLPOS$VIEW == 'periph' ,]

# mean error for each target per position x & y - not absolute
XPOS <- summarySEwithin(data = ALLPERIPH, measurevar = 'LANDx_mm', 
                              betweenvars = 'DIAGNOSIS', withinvars = c('ECCx', 'VIEW','SIDE'))
names(XPOS)[8] <- 'SEX'
YPOS <- summarySEwithin(data = ALLPERIPH, measurevar = 'LANDy_mm', 
                              betweenvars = 'DIAGNOSIS', withinvars = c('ECCx','VIEW','SIDE'))
names(YPOS)[8] <- 'SEY'

# combining
ALLPOS <- merge(XPOS, YPOS, 
                      by = c('DIAGNOSIS', 'VIEW', 'SIDE', 'N', 'ECCx'))
ALLPOS$DIAGNOSIS <- factor(ALLPOS$DIAGNOSIS, levels = c('HC','MCI','AD'))
ALLPOS$ECCy <- 0
ALLPOS$ECCx <- as.numeric(as.character(ALLPOS$ECCx))

for (l in 1:length(ALLPOS$ECCx)) {
  if (isTRUE(ALLPOS$SIDE[l] == 'left')){
    ALLPOS$ECCx[l] <- (ALLPOS$ECCx*(-1))[l]
  }
}

###### PLOTTING X BY Y ######
# peripheral -left
ggplot(ALLPOS[ALLPOS$SIDE=='left' ,]) +
  geom_point(aes(LANDx_mm, LANDy_mm, shape = DIAGNOSIS, colour = as.factor(ECCx)),
             size = 4.5) +
  geom_point(aes(as.numeric(ECCx), ECCy, colour= as.factor(ECCx)), shape = 21, size = 10) +
#  geom_errorbar(aes(x = LANDx_mm, ymin=LANDy_mm-ci.y, ymax=LANDy_mm+ci.y, 
#                    colour = DIAGNOSIS), width = .17, size = .75) +
  geom_errorbarh(aes(xmin=LANDx_mm-ci.x, xmax=LANDx_mm+ci.x, y = LANDy_mm,
                     colour = as.factor(ECCx)), height = .1, size = .75) +
  geom_hline(yintercept = 0, size = .75, linetype = 'dashed') + 
  scale_colour_manual(values = c('grey50','darkblue','darkgreen')) +
  ylim(-6,2) + xlim(-130,1) +
  labs(x = 'Reach endpoint (x-axis, mm)', y = 'Reach endpoint (y-axis, mm)', 
       title = 'Left') +
  theme_classic() + theme(legend.position = 'none',
                          legend.title = element_blank(),
                          legend.direction = 'horizontal') -> LPP

# peripheral -right
ggplot(ALLPOS[ALLPOS$SIDE=='right' ,]) +
  geom_point(aes(LANDx_mm, LANDy_mm, shape = DIAGNOSIS, colour = as.factor(ECCx)),
             size = 4.5) +
  geom_point(aes(as.numeric(ECCx), ECCy, colour= as.factor(ECCx)), shape = 21, size = 10) +
  #  geom_errorbar(aes(x = LANDx_mm, ymin=LANDy_mm-ci.y, ymax=LANDy_mm+ci.y, 
  #                    colour = DIAGNOSIS), width = .17, size = .75) +
  geom_errorbarh(aes(xmin=LANDx_mm-ci.x, xmax=LANDx_mm+ci.x, y = LANDy_mm,
                     colour = as.factor(ECCx)), height = .1, size = .75) +
  geom_hline(yintercept = 0, size = .75, linetype = 'dashed') + 
  scale_colour_manual(values = c('grey50','darkblue','darkgreen')) +
  ylim(-6,2) + xlim(-1,130) +
  labs(x = 'Reach endpoint (x-axis, mm)', y = '', 
       title = 'Right') +
  theme_classic() + theme(legend.position = c(.75,.95),
                          legend.title = element_blank(),
                          legend.direction = 'horizontal') -> RPP

## combining all
TWOD <- ggarrange(LPP,RPP,
                  ncol = 2, nrow = 1,
                  widths = c(1,1),
                  hjust = .1)
TWOD
# saving
ggsave('2DREACH.png', plot = last_plot(),  device = NULL, dpi = 300, 
       width = 8, height = 4, path = anaPath)

