# EXPLORATORY ANALYSIS FOR RADIAL REACHING TASK 
## AG.Mitchell 22.04.20 ##
library(readr)
library(ggplot2)
library(Rmisc)
library(gridExtra)
library(plyr)
library(tidyverse)
library(reshape2)
library(Hmisc)
library(ggpubr)
library(ez)
library(psychReport)

###### GETTING DATA ######
#on mac
anaPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/analysis/radial_reaching'
#dataPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/data'
#on pc
#dataPath <- 'S:/groups/DMT/data'
#anaPath <- 'S:/groups/DMT/analysis/radial_reaching'
setwd(anaPath)

res <- read.csv('all_radial-reaching_compiled.csv')
## find outliers and remove 
xclude <- read.csv('radial-reaching_outliers.csv')
res <- res[!(res$PPT %in% xclude$PPT), ]
# removing outer target
res <- res[res$POSITION != 400 & res$POSITION != -400 ,]

##### DIRECTIONAL ERROR: PMI #####
dir_medians <- aggregate(LANDx ~ PPT * VIEW * SIDE * POSITION * SITE * GRP * DIAGNOSIS, 
                         median, data = res)
# remove outlier in P101
dir_medians <- dir_medians[!(dir_medians$PPT == '101' & dir_medians$LANDx > 50) ,]

dir_means <- aggregate(LANDx ~ PPT* VIEW * SIDE * SITE * GRP * DIAGNOSIS, 
                       mean, data = dir_medians)

# casting by task
dPMIdata <- dcast(dir_means, PPT+DIAGNOSIS+GRP+SIDE+SITE ~ VIEW)
dPMIdata$PMI <- dPMIdata$Peripheral - dPMIdata$Free
write.csv(dPMIdata, 'radial-reaching_dirPMI.csv', row.names = FALSE)

# plotting this
ggplot(dir_means, aes(x = VIEW, y = LANDx, colour = DIAGNOSIS)) +
  geom_point(shape = 1, size = 2) +
  geom_line(aes(group = PPT), size = 0.5, alpha = .5) +
  facet_grid(cols = vars(SIDE), rows = vars(DIAGNOSIS)) + 
  labs(title = 'Radial reaching',
       x = '', y = 'Directional error (mm)', element_text(size = 12)) +
  theme_bw() + theme(legend.position = 'none', text = element_text(size = 10),
                     strip.text.x = element_text(size = 8)) 

# dPMI plot 
ggplot(dPMIdata, aes(x = SIDE, y = PMI, colour = DIAGNOSIS)) + 
  geom_point(shape = 1, size = 1.5, stroke = .8) +
  geom_line(aes(group = PPT), alpha = .5, size = .5) +
  stat_summary(aes(y = PMI, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1, size = 3, group = 1) +
  labs(title = 'Radial Reaching', x = 'Side', y = 'dPMI (mm)', 
       element_text(size = 12)) +
  facet_wrap(~DIAGNOSIS) +
  theme_bw() + theme(legend.position = 'none', text = element_text(size = 10),
                     strip.text.x = element_text(size = 10))

ggsave('radial-dirPMI.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

## ANOVA ##
dPMIanova <- dPMIdata[dPMIdata$PPT != 212 & 
                        dPMIdata$PPT != 310 ,] #removing participants where we only have 1 data-point

# FULL ANOVA ON FILTERED DATA
DPMI_ANOVA <- ezANOVA(
  data = dPMIanova
  , dv = .(PMI)
  , wid = .(PPT)
  , within = .(SIDE)
  , between = .(DIAGNOSIS)
  , type = 3,
  return_aov = TRUE,
  detailed = TRUE
)

DPMI_ANOVA$ANOVA
DPMI_ANOVA$`Mauchly's Test for Sphericity`
DPMI_ANOVA$`Sphericity Corrections`
aovDPMI <- aovEffectSize(ezObj = DPMI_ANOVA, effectSize = "pes")
aovDispTable(aovDPMI)

##### DIRECTIONAL ERROR: ECCENTRICITY #####
# adding eccentricity as a function of side to medians
dir_medians$ECC <- dir_medians$POSITION
#making left side negative
dir_medians$POSITION <- abs(dir_medians$ECC)
dir_medians$ECC <- factor(dir_medians$ECC)
dir_medians$POSITION <- factor(dir_medians$POSITION)

## plotting data frame
av_ecc <- summarySE(dir_medians, measurevar = 'LANDx', 
                    groupvar = c('DIAGNOSIS','ECC','VIEW','SIDE'), na.rm = TRUE)
av_ecc$DIAGNOSIS <- factor(av_ecc$DIAGNOSIS, levels = c('HC','MCI','AD'))
# plot 
ggplot(av_ecc, aes(x = ECC, y = LANDx, group = DIAGNOSIS, colour = DIAGNOSIS)) +
  geom_point(size = 3, position = position_dodge(width = .4)) +
  geom_errorbar(aes(ymin=LANDx-ci, ymax=LANDx+ci), 
                width=.4, position = position_dodge(width = .4)) + 
  geom_line(aes(group = DIAGNOSIS), size = 0.7, position = position_dodge(width = .4)) +
  geom_hline(yintercept = 0) + 
  labs(title= 'Radial reaching', x = 'Eccentricity (mm)', y = 'Directional error (mm)') +
  facet_grid(~VIEW) + theme_classic() +
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12)
  )

ggsave('Radial-dAE-ECC.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

## ANOVA ##
# remove participants with unbalanced data
dECCanova <- dir_medians[dir_medians$PPT != 212 & 
                           dir_medians$PPT != 310 &
                           dir_medians$PPT != 101
                         ,]

# FULL ANOVA ON ECCENTRICITY DATA
DECC_ANOVA <- ezANOVA(
  data = dECCanova
  , dv = .(LANDx)
  , wid = .(PPT)
  , within = .(VIEW, SIDE, POSITION)
  , between = .(DIAGNOSIS)
  , type = 3,
  return_aov = TRUE,
  detailed = TRUE
)

DECC_ANOVA$ANOVA
DECC_ANOVA$`Mauchly's Test for Sphericity`
DECC_ANOVA$`Sphericity Corrections`
aovDECC <- aovEffectSize(ezObj = DECC_ANOVA, effectSize = "pes")
aovDispTable(aovDECC)

###### MOVEMENT TIME #######
MT_medians <- aggregate(MT ~ PPT * VIEW * SIDE * DOM * POSITION * SITE * GRP * DIAGNOSIS, 
                        median, data = res)
MT_means <- aggregate(MT ~ PPT * VIEW * SIDE * DOM * SITE * GRP * DIAGNOSIS,
                      mean, data = MT_medians)
# average across side
MTav <- aggregate(MT ~ PPT * VIEW * SITE * GRP * DIAGNOSIS,
                  mean, data = MT_medians)
## plotting :)
# both sides
ggplot(MT_means, aes(x = DOM, y = MT, colour = SITE, group = DIAGNOSIS)) +
  geom_point(shape = 1, size = 2) +
  geom_line(aes(group = PPT), size = 0.5, alpha = .5) +
  facet_grid(cols = vars(VIEW), rows = vars(DIAGNOSIS)) + ylim(0, 1100) +
  labs(title = 'Radial reaching', x = 'Side', 
       y = 'Reach duration (ms)', element_text(size = 12)) +
  theme_bw() + theme(legend.position = 'none', 
                     text = element_text(size = 10),
                     strip.text.x = element_text(size = 10)) 

# average across sides
MTav$DIAGNOSIS <- factor(MTav$DIAGNOSIS, levels = c('HC','MCI','AD'))
ggplot(MTav, aes(x = VIEW, y = MT, colour = SITE, group = PPT)) +
  geom_point(shape = 16, size = 2, position = position_dodge(width = .3)) +
  geom_line(aes(group = PPT), size = 0.5, alpha = .5, 
            position = position_dodge(width = .3)) +
  stat_summary(aes(y = MT, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1, size = 4, group = 1) +
  facet_wrap(~DIAGNOSIS) + 
  labs(title = 'Radial reaching', x = '', 
       y = 'Reach duration (ms)', element_text(size = 12)) +
  theme_classic() + theme(legend.position = 'bottom', 
                          text = element_text(size = 10),
                          strip.text.x = element_text(size = 10)
  ) 

ggsave('radial-MT.png', plot = last_plot(),  device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

# by eccentricity
MT_medians$ECC <- MT_medians$POSITION
#making left side negative
MT_medians$POSITION <- abs(MT_medians$ECC)
MT_medians$ECC <- factor(MT_medians$ECC)
MT_medians$POSITION <- factor(MT_medians$POSITION)

# summary data
MTecc <- summarySE(MT_medians, measurevar = 'MT', 
                   groupvar = c('DIAGNOSIS','ECC','VIEW','SIDE'), na.rm = TRUE)
MTecc$DIAGNOSIS <- factor(MTecc$DIAGNOSIS, levels = c('HC','MCI','AD'))

ggplot(MTecc, aes(x = ECC, y = MT, group = DIAGNOSIS, colour = VIEW)) +
  geom_point(size = 3, position = position_dodge(width = .4)) +
  geom_errorbar(aes(ymin=MT-ci, ymax=MT+ci), 
                width=.4, position = position_dodge(width = .4)) + 
  geom_line(aes(group = VIEW), size = 0.7, position = position_dodge(width = .4)) +
  labs(title = 'Radial reaching',
       x = 'Eccentricity (°)', y = 'Movement time (ms)') +
  facet_grid(~DIAGNOSIS) + theme_classic() +
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12)
  )

ggsave('radial-MTecc.png', plot = last_plot(),  device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

## ANOVA ## 
MTanova <- MT_means[MT_means$PPT != 212 & MT_means$PPT != 310 ,]

# FULL ANOVA ON MOVEMENT TIME DATA
MT_ANOVA <- ezANOVA(
  data = MTanova
  , dv = .(MT)
  , wid = .(PPT)
  , within = .(VIEW, DOM)
  , between = .(DIAGNOSIS)
  , between_covariates = .(SITE)
  , type = 3,
  return_aov = TRUE,
  detailed = TRUE
)

MT_ANOVA$ANOVA
MT_ANOVA$`Mauchly's Test for Sphericity`
MT_ANOVA$`Sphericity Corrections`
aovMT <- aovEffectSize(ezObj = MT_ANOVA, effectSize = "pes")
aovDispTable(aovMT)

# FULL ANOVA ON MOVEMENT TIME DATA BY ECCENTRICITY
MTECCanova <- MT_medians[MT_medians$PPT != 212 & 
                           MT_medians$PPT != 310 &
                           MT_medians$PPT != 315 &
                           MT_medians$PPT != 302 ,]

MTECCanova$POSITION <- factor(MTECCanova$POSITION)
MTECC_ANOVA <- ezANOVA(
  data = MTECCanova
  , dv = .(MT)
  , wid = .(PPT)
  , within = .(VIEW, DOM, POSITION)
  , between = .(DIAGNOSIS)
  , between_covariates = .(SITE)
  , type = 3,
  return_aov = TRUE,
  detailed = TRUE
)

MTECC_ANOVA$ANOVA
MTECC_ANOVA$`Mauchly's Test for Sphericity`
MTECC_ANOVA$`Sphericity Corrections`
aovMTECC <- aovEffectSize(ezObj = MTECC_ANOVA, effectSize = "pes")
aovDispTable(aovMTECC)



# all target positions - by side
ggplot(res_rt_means, aes(x = VIEW, y = RT, colour = GRP), position = position_dodge(.2)) + 
  geom_point(shape = 1, size = 1.5, stroke = .8) +
  facet_grid(cols = vars(SIDE), rows = vars(GRP)) +
  geom_line(aes(group = PPT), alpha = .5, size = .5) +
  scale_colour_manual(values = c('black', 'grey40')) +
  ylim(0,1000) + labs(title = 'Radial reaching', x = 'View', y = 'Reaction time (ms)', 
                        element_text(size = 12)) +
  theme_bw() + theme(legend.position = 'none', text = element_text(size = 10),
                     strip.text.x = element_text(size = 10)) -> RTplot

ggsave('RT_side.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)


# plot for everything
res_rt_meansall <- aggregate(RT~VIEW*PPT*GRP, mean, data = res_rt_means)

ggplot(res_rt_meansall, aes(x = VIEW, y = RT, colour = GRP), position = position_dodge(.2)) + 
  geom_point(shape = 1, size = 1.5, stroke = .8) +
  facet_wrap(~GRP) +
  geom_line(aes(group = PPT), alpha = .5, size = .5) +
  scale_colour_manual(values = c('black', 'grey40')) +
  ylim(0,1000) + labs(title = 'Radial reaching', x = 'Viewing condition', y = 'Reaction time (ms)', 
                        element_text(size = 12)) +
  theme_bw() + theme(legend.position = 'none', text = element_text(size = 10),
                     strip.text.x = element_text(size = 10)) -> RTplot

ggsave('RT.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)



###### normalised movement time after peak speed
res$NMTPS <- (res$MT - res$TPS)/res$MT
plotNMTPS <- summarySE(data=res, measurevar = "NMTPS", 
                       groupvars = c("GRP", "POSITION", "VIEW"), na.rm = TRUE)
plotNMTPS <- na.omit(plotNMTPS)

# plot
ggplot(plotNMTPS, aes(x=POSITION, y=NMTPS, colour=GRP, group=GRP)) +
  geom_point(size=5, alpha=.5, position=position_dodge(width=.3)) +
  geom_errorbar(aes(ymin=NMTPS-ci, ymax=NMTPS+ci), width=.4, position=position_dodge(width=.3)) +
  geom_line(position=position_dodge(width=.3)) +
  facet_wrap(~VIEW) +
  theme_bw()

ggsave('normMTafterPS.png', plot = last_plot(), device = NULL, dpi = 300, 
       width = 8, height = 7, path = anaPath)

### correlate reaching error with NMTPS
# creating relevant data-frame
NMTPS_mean <- aggregate(NMTPS ~  VIEW*SIDE*PPT*GRP, mean, data = res)
NMTPS_mean <- dcast(NMTPS_mean, PPT+GRP+SIDE ~ VIEW)
# rename columns
levels(NMTPS_mean$SIDE) <- c('Left', 'Right')
levels(NMTPS_mean$GRP) <- c('Control', 'Patient') 
names(NMTPS_mean)[4] <- 'NMTPS_Free'
names(NMTPS_mean)[5] <- 'NMTPS_Periph'
NMTPS_mean$PPT <- substr(NMTPS_mean$PPT, 4, 6)
NMTPS_mean$COST <- NMTPS_mean$NMTPS_FREE - NMTPS_mean$NMTPS_PERIPH

# merge with PMI
test <- merge(PMIdata, NMTPS_mean, by = c('PPT','GRP', 'SIDE'), all = TRUE)
Periphav <- aggregate(Peripheral ~ PPT * GRP, mean, data = test)
NMTPSav <- aggregate(NMTPS_Free ~ PPT * GRP, mean, data = test)
testav <- merge(Periphav, NMTPSav, by = c('PPT', 'GRP'), all = TRUE)

# plot correlation
ggscatter(test, x = 'Peripheral', y = 'NMTPS_Free', add = 'reg.line', conf.int = TRUE,
          cor.coef = TRUE, cor.method = 'pearson') + 
  facet_grid(cols = vars(SIDE), rows = vars(GRP)) + 
  ylab('Normalised MTAPS (free)') + xlab('Peripheral reaching error (deg)')

# plot correlation - averaged across sides
ggscatter(testav, x = 'Peripheral', y = 'NMTPS_Free', add = 'reg.line', conf.int = TRUE,
          cor.coef = TRUE, cor.method = 'spearman') + 
  facet_wrap(~GRP) + 
  ylab('Normalised MTAPS (free)') + xlab('Peripheral reaching error (deg)')

ggsave('normMTAPS-accuracy_correlation.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

#### time to peak speed
plotTPS <- summarySE(data=res, measurevar = "TPS", 
                     groupvars = c("GRP", "POSITION", "VIEW"), na.rm = TRUE)
plotTPS <- na.omit(plotTPS)

ggplot(plotTPS, aes(x=POSITION, y=TPS, colour=GRP, group=GRP)) +
  geom_point(size=5, alpha=.5, position=position_dodge(width=.3)) +
  geom_errorbar(aes(ymin=TPS-ci, ymax=TPS+ci), width=.4, position=position_dodge(width=.3)) +
  geom_line(position=position_dodge(width=.3)) +
  facet_wrap(~VIEW) +
  theme_bw()

