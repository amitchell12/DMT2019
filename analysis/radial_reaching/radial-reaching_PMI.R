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

# summary data
res_medians <- aggregate(AEdeg ~ POSITION*VIEW*SIDE*PPT*GRP, mean, data = res)
colnames(res_medians)[colnames(res_medians)=='AEdeg'] <- 'AEmed'
res_means <- aggregate(AEmed ~ VIEW*SIDE*PPT*GRP, mean, data = res_medians)
colnames(res_means)[colnames(res_means)=='AEmed'] <- 'AEmean'

# changing levels of PMI for plotting
res_means$VIEW <- factor(res_means$VIEW) #changing so only 2 levels recorded
res_means$SIDE <- factor(res_means$SIDE, levels = c('LEFT', 'RIGHT'))
levels(res_means$SIDE) <- c('Left', 'Right')
levels(res_means$GRP) <- c('Control', 'Patient') #changing group name from 1 = control, 2 = AD
levels(res_means$VIEW) <- c('Free', 'Peripheral')
res_means$PPT <- substr(res_means$PPT, 4, 6)

# casting by task
PMIdata <- dcast(res_means, PPT+GRP+SIDE ~ VIEW)
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

#isolating peripheral data
res_periph <- res_means[res_means$VIEW == 'Peripheral' ,]

# PMI plot
ggplot(res_periph, aes(x = SIDE, y = AEmean, colour = GRP), position = position_dodge(.2)) + 
  geom_point(shape = 1, size = 4) +
  geom_line(aes(group = PPT), alpha = .5, size = .8) +
  scale_colour_manual(values = c('grey40', 'grey40')) +
  stat_summary(aes(y = AEmean, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1, size = 5, group = 1) +
  ylim(-.5,8) + labs(title = 'Radial reaching', x = 'Side', y = 'PMI (deg)', 
                     element_text(size = 14)) +
  facet_wrap(~GRP) +
  theme_bw() + theme(legend.position = 'none', text = element_text(size = 14),
                     strip.text.x = element_text(size = 12)) -> PMIplot

ggsave('periphmeans_plot.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, width = 7, height = 4, path = anaPath)

## summary data
meanPMI_side <- summarySE(PMIdata, measurevar = 'PMI', groupvar = c('GRP', 'SIDE'),
                          na.rm = TRUE)
meanPMI_all <- summarySE(PMIdata, measurevar = 'PMI', groupvar = c('GRP'),
                         na.rm = TRUE)

##### directional error calc #####
dir_medians <- aggregate(LANDx_deg ~ POSITION*VIEW*SIDE*PPT*GRP, mean, data = res)
colnames(dir_medians)[colnames(dir_medians)=='LANDx_deg'] <- 'dirmed'
dir_means <- aggregate(dirmed ~ VIEW*SIDE*PPT*GRP, mean, data = dir_medians)
colnames(dir_means)[colnames(dir_means)=='dirmed'] <- 'dirmean'

# changing levels of PMI for plotting
dir_means$VIEW <- factor(dir_means$VIEW) #changing so only 2 levels recorded
dir_means$SIDE <- factor(dir_means$SIDE, levels = c('LEFT', 'RIGHT'))
levels(dir_means$SIDE) <- c('Left', 'Right')
levels(dir_means$GRP) <- c('Control', 'AD') #changing group name from 1 = control, 2 = AD
levels(dir_means$VIEW) <- c('Free', 'Peripheral')
dir_means$PPT <- substr(dir_means$PPT, 4, 6)

# casting by task
dPMIdata <- dcast(dir_means, PPT+GRP+SIDE ~ VIEW)
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

