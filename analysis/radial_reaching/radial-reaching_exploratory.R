###### EXPLORATORY ANALYSIS FOR RADIAL REACHING TASK ######
## AG.Mitchell 22.04.20 ##

library(readr)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(Rmisc)

#on mac
#anaPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/analysis/radial_reaching'
#dataPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/data'
#on pc
dataPath <- 'S:/groups/DMT/data'
anaPath <- 'S:/groups/DMT/analysis/radial_reaching'
setwd(anaPath)

res <- read.csv('radial_reaching_compiled.csv')


##### summary response time #####
#response time
res_rt_posmean <- aggregate(RT ~ POSITION*VIEW*SIDE*PPT*GRP, mean, data = res)
res_rt_means <- aggregate(RT ~ VIEW*SIDE*PPT*GRP, mean, data = res)

res_rt_means$VIEW <- factor(res_rt_means$VIEW) #changing so only 2 levels recorded
res_rt_means$SIDE <- factor(res_rt_means$SIDE, levels = c('LEFT', 'RIGHT'))
levels(res_rt_means$SIDE) <- c('Left', 'Right')
levels(res_rt_means$GRP) <- c('Control', 'Patient') #changing group name from 1 = control, 2 = AD
levels(res_rt_means$VIEW) <- c('Free', 'Peripheral')
res_rt_means$PPT <- substr(res_rt_means$PPT, 4, 6)

#movement time
res_mt_posmean <- aggregate(MT ~ POSITION*VIEW*SIDE*PPT*GRP, mean, data = res)
res_mt_means <- aggregate(MT ~ VIEW*SIDE*PPT*GRP, mean, data = res)

res_mt_means$VIEW <- factor(res_mt_means$VIEW) #changing so only 2 levels recorded
res_mt_means$SIDE <- factor(res_mt_means$SIDE, levels = c('LEFT', 'RIGHT'))
levels(res_mt_means$SIDE) <- c('Left', 'Right')
levels(res_mt_means$GRP) <- c('Control', 'Patient') #changing group name from 1 = control, 2 = AD
levels(res_mt_means$VIEW) <- c('Free', 'Peripheral')
res_mt_means$PPT <- substr(res_mt_means$PPT, 4, 6)

# plotting movement time
# per target position (eccentricity)
ggplot(res_mt_posmean, aes(x = POSITION, y = MT, colour = GRP), position = position_dodge(.2)) + 
  geom_point(shape = 1, size = 1.5, stroke = .8) +
  facet_grid(cols = vars(VIEW), rows = vars(GRP)) +
  scale_colour_manual(values = c('grey40', 'grey40')) +
  stat_summary(aes(y = MT, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1, size = 2, group = 1) +
  ylim(400,1100) + labs(title = 'Radial reaching', x = 'Target position (mm)', y = 'Movement time (ms)', 
                    element_text(size = 12)) +
  theme_bw() + theme(legend.position = 'none', text = element_text(size = 10),
                     strip.text.x = element_text(size = 10)) -> MTPosplot

ggsave('MT_position.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

# all target positions - by side
ggplot(res_mt_means, aes(x = VIEW, y = MT, colour = GRP), position = position_dodge(.2)) + 
  geom_point(shape = 1, size = 1.5, stroke = .8) +
  facet_grid(cols = vars(SIDE), rows = vars(GRP)) +
  geom_line(aes(group = PPT), alpha = .5, size = .5) +
  scale_colour_manual(values = c('black', 'grey40')) +
  ylim(400,1100) + labs(title = 'Radial reaching', x = 'View', y = 'Movement time (ms)', 
                        element_text(size = 12)) +
  theme_bw() + theme(legend.position = 'none', text = element_text(size = 10),
                     strip.text.x = element_text(size = 10)) -> MTplot

ggsave('MT_side.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)


# plot for everything
res_mt_meansall <- aggregate(MT~VIEW*PPT*GRP, mean, data = res_mt_means)

ggplot(res_mt_meansall, aes(x = VIEW, y = MT, colour = GRP), position = position_dodge(.2)) + 
  geom_point(shape = 1, size = 1.5, stroke = .8) +
  facet_wrap(~GRP) +
  geom_line(aes(group = PPT), alpha = .5, size = .5) +
  scale_colour_manual(values = c('black', 'grey40')) +
  ylim(400,1000) + labs(title = 'Radial reaching', x = 'Viewing condition', y = 'Movement time (ms)', 
                        element_text(size = 12)) +
  theme_bw() + theme(legend.position = 'none', text = element_text(size = 10),
                     strip.text.x = element_text(size = 10)) -> MTplot

ggsave('MT.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

# plotting reaction time
# per target position (eccentricity)
ggplot(res_rt_posmean, aes(x = POSITION, y = RT, colour = GRP), position = position_dodge(.2)) + 
  geom_point(shape = 1, size = 1.5, stroke = .8) +
  facet_grid(cols = vars(VIEW), rows = vars(GRP)) +
  scale_colour_manual(values = c('grey40', 'grey40')) +
  stat_summary(aes(y = RT, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1, size = 2, group = 1) +
  ylim(0,1000) + labs(title = 'Radial reaching', x = 'Target position (mm)', y = 'Reaction time (ms)', 
                        element_text(size = 12)) +
  theme_bw() + theme(legend.position = 'none', text = element_text(size = 10),
                     strip.text.x = element_text(size = 10)) -> RTPosplot

ggsave('RT_position.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

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

