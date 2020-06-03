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

res <- read.csv('radial-reaching_compiled.csv')
# adding site as a factor 
res$PPT <- substr(res$PPT, 4, 6)
res$SITE <- factor(substr(res$PPT, 1, 1)) # 1/2 edinburgh, 3/4 norwich
for (i in 1:length(res$SITE)){
  if (isTRUE(res$SITE[i] == "2")){
    res$SITE[i] = 1
  }
}
# changing levelsfor plotting
res$VIEW <- factor(res_means$VIEW) #changing so only 2 levels recorded
res$SIDE <- factor(res_means$SIDE, levels = c('LEFT', 'RIGHT'))
levels(ress$SIDE) <- c('Left', 'Right')
res$GRP <- factor(res$GRP)
levels(res$GRP) <- c('Control', 'Patient') #changing group name from 1 = control, 2 = AD
levels(res$VIEW) <- c('Free', 'Peripheral')

# summary data
res_medians <- aggregate(AEdeg ~ PPT*POSITION*VIEW*SIDE*GRP*SITE, mean, data = res)
colnames(res_medians)[colnames(res_medians)=='AEdeg'] <- 'AEmed'
res_means <- aggregate(AEmed ~ PPT*VIEW*SIDE*GRP*SITE, mean, data = res_medians)
colnames(res_means)[colnames(res_means)=='AEmed'] <- 'AEmean'

# casting by task
PMIdata <- dcast(res_means, PPT+GRP+SIDE+SITE ~ VIEW)
PMIdata$PMI <- PMIdata$Peripheral - PMIdata$Free
write.csv(PMIdata, 'radial-reaching_PMI.csv', row.names = FALSE)


##### plots #####
# mean plot 
ggplot(res_means, aes(x = SIDE, y = AEmean, colour = GRP)) +
  geom_point(shape = 1, size = 2) +
  geom_line(aes(group = PPT), size = 0.5, alpha = .5) +
  facet_grid(cols = vars(VIEW), rows = vars(GRP)) + ylim(-0.5,7) +
  labs(x = 'Side', y = 'Mean AE (deg)', element_text(size = 12)) +
  scale_colour_manual(values = c('black', 'grey50')) +
  theme_bw() + theme(legend.position = 'none', text = element_text(size = 10),
                     strip.text.x = element_text(size = 8)) -> meansPlot

ggsave('allmeans_plot.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

# PMI plot
ggplot(PMIdata, aes(x = SIDE, y = PMI, colour = GRP), position = position_dodge(.2)) + 
  geom_point(shape = 1, size = 4) +
  geom_line(aes(group = PPT), alpha = .5, size = .8) +
  scale_colour_manual(values = c('grey40', 'grey40')) +
  stat_summary(aes(y = PMI, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1, size = 5, group = 1) +
  ylim(-.5,8) + labs(title = 'Radial reaching', x = 'Side', y = 'PMI (deg)', 
                     element_text(size = 14)) +
  facet_wrap(~GRP) +
  theme_bw() + theme(legend.position = 'none', text = element_text(size = 14),
                     strip.text.x = element_text(size = 12)) -> PMIplot

ggsave('PMIplot.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, width = 7, height = 4, path = anaPath)

## summary data
meanPMI_side <- summarySE(PMIdata, measurevar = 'PMI', groupvar = c('GRP', 'SIDE'),
                          na.rm = TRUE)
meanPMI_all <- summarySE(PMIdata, measurevar = 'PMI', groupvar = c('GRP'),
                         na.rm = TRUE)

##### directional error calc #####
dir_medians <- aggregate(LANDx_deg ~ PPT*POSITION*VIEW*SIDE*GRP*SITE, mean, data = res)
colnames(dir_medians)[colnames(dir_medians)=='LANDx_deg'] <- 'dirmed'
dir_means <- aggregate(dirmed ~ PPT*VIEW*SIDE*GRP*SITE, mean, data = dir_medians)
colnames(dir_means)[colnames(dir_means)=='dirmed'] <- 'dirmean'

# casting by task
dPMIdata <- dcast(dir_means, PPT+GRP+SIDE+SITE ~ VIEW)
dPMIdata$PMI <- dPMIdata$Peripheral - PMIdata$Free
write.csv(dPMIdata, 'radial-reaching_dPMI.csv', row.names = FALSE)

## plotting
# mean plot 
ggplot(dir_means, aes(x = VIEW, y = dirmean, colour = GRP)) +
  geom_point(shape = 1, size = 2) +
  geom_line(aes(group = PPT), size = 0.5, alpha = .5) +
  facet_grid(cols = vars(SIDE), rows = vars(GRP)) + ylim(-8,8) +
  labs(x = 'Side', y = 'Directional error (deg)', element_text(size = 12)) +
  scale_colour_manual(values = c('black', 'grey50')) +
  theme_bw() + theme(legend.position = 'none', text = element_text(size = 10),
                     strip.text.x = element_text(size = 8)) -> meansPlot

ggsave('directionalerror_plot.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

# PMI plot
ggplot(dPMIdata, aes(x = SIDE, y = PMI, colour = GRP), position = position_dodge(.2)) + 
  geom_point(shape = 1, size = 1.5, stroke = .8) +
  geom_line(aes(group = PPT), alpha = .5, size = .5) +
  scale_colour_manual(values = c('grey40', 'grey40')) +
  stat_summary(aes(y = PMI, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1, size = 2, group = 1) +
  ylim(-8,8) + labs(title = 'Radial reaching', x = 'Side', y = 'dPMI (deg)', 
                     element_text(size = 12)) +
  facet_wrap(~GRP) +
  theme_bw() + theme(legend.position = 'none', text = element_text(size = 10),
                     strip.text.x = element_text(size = 8)) -> PMIplot

ggsave('dPMI_plot.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

##### outlier calculation, for controls only ######
controlData <- PMIdata[PMIdata$PPT < 200, ]

# median values for each side
tmp <- aggregate(controlData$PMI,  by=list(SIDE = controlData$SIDE), FUN=median)
names(tmp)[2] <- 'med'
controlData <- merge(tmp, controlData)

# calculating MAD for each (absolute value)
controlData$AD <- abs(controlData$PMI - controlData$med)
tmp <- aggregate(controlData$AD,  by=list(side = controlData$SIDE), FUN=median)
names(tmp)[2] <- 'MAD'
controlData <- merge(controlData, tmp)

# adjusted z-score from these values
controlData$az <- (controlData$PMI - controlData$med)/(controlData$MAD * 1.4826)
controlData$z <- scale(controlData$PMI)
controlData$PPT <- factor(controlData$PPT)

plot_name = 'adjustedZ.png'
AZplot <- ggplot(controlData, aes(x = SIDE, y = az, colour = PPT)) +
  geom_point(size = 3, position = position_dodge(.1)) + ylim(-3,6) + 
  stat_summary(aes(y = az, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', group = 1) + 
  labs(x = '', y = 'Adjusted z-score', element_text(size = 13)) +
  theme_bw() + theme(legend.position = "right", legend.title = element_blank())

ggsave(plot_name, plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, width = 5, height = 6.5, path = anaPath)

####### outlier removal #######
# find controls with az > 2.5 and flag-up in PMI data-set
for (l in 1:length(controlData$PPT)){
  if (isTRUE(controlData$az[l] > '2.5')){
    controlData$outlier[l] = 1
  }
  else 
    controlData$outlier[l] = 0
}
#reorder control data
controlData <- controlData[order(controlData$PPT), ]

# adding to PMI data set, then removing
PMIfilter <- merge(PMIdata, controlData, all = TRUE)
PMIfilter <- PMIfilter[, c(1:7,14)]
# save with outlier information
write.csv(PMIfilter, 'radial-outliers.csv', row.names = FALSE)

## removing outliers
for (l in 49:length(PMIfilter$outlier)){
  PMIfilter$outlier[l] = 0
}
PMIfilter <- PMIfilter[PMIfilter$outlier < 1, ]

### PLOT FILTERED PMI DATA
ggplot(PMIfilter, aes(x = SIDE, y = PMI, colour = GRP), position = position_dodge(.2)) + 
  geom_point(shape = 1, size = 4) +
  geom_line(aes(group = PPT), alpha = .5, size = .8) +
  scale_colour_manual(values = c('grey50', 'black')) +
  stat_summary(aes(y = PMI, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1, size = 5, group = 1) +
  ylim(-.5,10) + labs(title = 'Radial Reaching', x = 'Side', y = 'Reaching error (deg)', 
                      element_text(size = 14)) +
  facet_wrap(~GRP) +
  theme_bw() + theme(legend.position = 'none', text = element_text(size = 14),
                     strip.text.x = element_text(size = 12)) -> PMIf_plot

ggsave('radialPMI-filtered.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, width = 7, height = 4, path = anaPath)

meanFPMI <- summarySE(PMIfilter, measurevar = 'PMI', groupvar = c('GRP', 'SIDE'),
                      na.rm = TRUE)
meanFPMI_all <- summarySE(PMIfilter, measurevar = 'PMI', groupvar = c('GRP'),
                          na.rm = TRUE)

#average across side
PMIfilter_av <- aggregate(PMI ~ PPT * SITE * GRP, mean, data = PMIfilter)
jitter <- position_jitter(width = 0.1, height = 0.1)

ggplot(PMIfilter_av, aes(x = GRP, y = PMI)) + 
  geom_point(position = jitter, shape = 21, size = 3) +
  stat_summary(aes(y = PMI, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1, size = 4, group = 1) +
  ylim(-.5,7) + labs(title = '', x = '', y = 'Reaching error (deg)', 
                     element_text(size = 8)) +
  theme_bw() + theme(legend.position = 'none', text = element_text(size = 10),
                     strip.text.x = element_text(size = 8)) -> PMIf_plot

ggsave('radialPMI-filtered-av.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, width = 3, height = 3, path = anaPath)


########### next steps: comparing patients to controls
