## A.G.Mitchell 26.05.2021
## SUPPLEMENTARY RAD ANALYSIS - X,Y coordinates at peak speed
## is there a difference between groups?

library(ggplot2)
library(ez)
library(psychReport)
library(Rmisc)

# loading files
anaPath <- '/Users/alex/Library/Mobile Documents/com~apple~CloudDocs/Documents/DMT/analysis/radial_reaching'
dataPath <- '/Users/alex/Library/Mobile Documents/com~apple~CloudDocs/Documents/DMT/data'

res <- read.csv('RAD_PSposition.csv')

# plotting to be sure
# plotting to get a look at data
res$POSITION <- factor(res$POSITION)
ggplot(res) + geom_point(aes(x = calx, y = caly, colour = POSITION), shape = 3) +
  geom_point(aes(x = px, y = py, colour = POSITION)) +
  facet_wrap(~PPT*VIEW)
ggsave('RAD_PSerr.png', plot = last_plot(), device = NULL, 
       path = anaPath, scale = 1, width = 15, height = 10, units = 'in')

###### CALCULATING VARS ######
# subtracting cal from reach endpoint
res$LANDx <- res$px - res$calx
res$LANDy <- res$py - res$caly

# absolute error - in mm
res$AE <- sqrt(res$LANDx^2 + res$LANDy^2) #mm

## angular error and amplitude error
# calculating target and end-point position relative to start-point
res$rX <- res$px - res$sX
res$rY <- res$py - res$sY
res$tX <- res$calx - res$sX
res$tY <- res$caly - res$sY

#recode response as target-relative ERRORS in polar coordinates
res$tANG <- (atan(res$tX/res$tY))*(180/pi)
res$tAMP <- sqrt(res$tX^2 + res$tY^2)
res$rANG <- (atan(res$rX/res$rY))*(180/pi)
res$rAMP <- sqrt(res$rX^2 + res$rY^2)

## calculating angular error and amplitude error
res$ANG_ERR <- res$rANG - res$tANG
res$AMP_ERR <- res$rAMP - res$tAMP

# plotting angular error
res$PPT <- factor(res$PPT)
ggplot(res) +
  geom_point(aes(x = POSITION, y = ANG_ERR, colour = PPT)) +
  facet_wrap(~DIAGNOSIS*VIEW)

###### ANG ERR ANOVA ######
res$ECC <- factor(abs(as.numeric(as.character(res$POSITION)))) #adding eccentricity = absolute target positio
res_medians <- aggregate(ANG_ERR ~ PPT * VIEW * SIDE * SITE * GRP * DIAGNOSIS * AGE * POSITION, 
                         median, data = res)
resANGERR <- aggregate(ANG_ERR ~ PPT*POSITION*VIEW*DIAGNOSIS*SITE*AGE, 
                   mean, data = res_medians)

resANGERR$POSITION <- factor(resANGERR$POSITION)
resANGERR$DIAGNOSIS <- factor(resANGERR$DIAGNOSIS)
resANGERR$AGE <- factor(resANGERR$AGE)
resANGERR$VIEW <- factor(resANGERR$VIEW)
resANGERR$SITE <- factor(resANGERR$SITE)

resANGERR <- resANGERR[!resANGERR$PPT == 212 & 
                  !resANGERR$PPT == 315 &
                  !resANGERR$PPT == 403 &
                  !resANGERR$PPT == 407,]

## CALCULATING ANOVA ##
ANGERR_ANOVA <- ezANOVA(
  data = resANGERR
  , dv = .(ANG_ERR)
  , wid = .(PPT)
  , within = .(VIEW, POSITION)
  , between = .(DIAGNOSIS)
  , between_covariates = .(SITE, AGE)
  , type = 3,
  return_aov = TRUE,
  detailed = TRUE
)

ANGERR_ANOVA$ANOVA
ANGERR_ANOVA$`Mauchly's Test for Sphericity`
ANGERR_ANOVA$`Sphericity Corrections`
aovANGERR <- aovEffectSize(ANGERR_ANOVA, effectSize = "pes")
aovDispTable(aovANGERR)

###### PLOTTING ######
# By eccentricity
ANGsummary <- summarySE(resANGERR, measurevar = 'ANG_ERR', 
                        groupvar = c('DIAGNOSIS', 'VIEW', 'POSITION'), na.rm = TRUE)
ANGsummary$DIAGNOSIS <- factor(ANGsummary$DIAGNOSIS, levels = c('HC','MCI','AD'))

ggplot(ANGsummary, aes(x = POSITION, y = ANG_ERR, group = DIAGNOSIS, colour = DIAGNOSIS, 
                       shape = DIAGNOSIS)) +
  geom_point(size = 3, position = position_dodge(.4)) +
  geom_errorbar(aes(ymin=ANG_ERR-ci, ymax=ANG_ERR+ci), 
                width=.4, position = position_dodge(.4)) + 
  geom_line(aes(group = DIAGNOSIS), size = 0.7, position = position_dodge(.4)) +
  scale_color_manual(values = c('black','grey30','grey60')) +
  labs(x = 'Eccentricity (mm)', y = 'Angular error at PS (°)') +
  facet_grid(~VIEW) + theme_classic() +
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12)
  ) -> ANGecc
ANGecc
ggsave('RAD_PSerr.png', plot = last_plot(),  device = NULL, dpi = 300, 
       width = 8, height = 6, path = anaPath)

# Mean for each participant
ANGERR_ECC <- aggregate(ANG_ERR ~ PPT * VIEW * POSITION * SIDE * SITE * DIAGNOSIS *AGE, 
                    median, data = res_medians)
ANGERR_means <- aggregate(ANG_ERR ~ PPT * VIEW * SITE * SIDE * DIAGNOSIS * AGE,
                      mean, data = ANGERR_ECC)
ANGERRgrp_means <- summarySE(ANGERR_means, measurevar = 'ANG_ERR', 
                             groupvars = c('DIAGNOSIS','VIEW','SIDE'))

ANGERR_means$DIAGNOSIS <- factor(ANGERR_means$DIAGNOSIS, levels = c('HC','MCI','AD'))

ggplot(ANGERR_means, aes(x = VIEW, y = ANG_ERR, colour = SITE, group = PPT)) +
  geom_point(shape = 16, size = 2, position = position_dodge(width = .3)) +
  geom_line(aes(group = PPT), size = 0.5, alpha = .5, 
            position = position_dodge(width = .3)) +
  scale_color_manual(values = c('dodgerblue3','grey50')) +
  stat_summary(aes(y = ANG_ERR, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1, size = 4, group = 1) +
    facet_wrap(~SIDE*DIAGNOSIS) + 
  labs(x = '', y = 'Angular error at PS (°)', element_text(size = 12)) +
  theme_classic() + theme(legend.position = 'bottom', 
                          axis.text = element_text(size = 10),
                          axis.title = element_text(size = 12),
                          strip.text = element_text(size = 12)
  ) 
ggsave('RAD_PSerr_means.png', plot = last_plot(),  device = NULL, dpi = 300, 
width = 6, height = 6, path = anaPath)


#pair-wise t-test
ANGERR_ECC$DIAGNOSIS <- factor(ANGERR_ECC$DIAGNOSIS)
ANGttest <- pairwise.t.test(ANGERR_ECC$ANG_ERR, ANGERR_ECC$DIAGNOSIS, p.adj = 'bonf')
print(ANGttest)
