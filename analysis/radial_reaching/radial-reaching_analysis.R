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
setwd(dataPath)

files <- list.files(path=dataPath, pattern = "*.TRJ", full.names = TRUE, recursive = TRUE)
idfiles <- list.files(path=dataPath, pattern = "*.dmt", full.names = TRUE, recursive = TRUE)

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

##### variables #####
visAngle <- function(size, distance){
  # this function calculates visual angle
  # size and distance must be in the same units
  Rad = 2*atan(size/(2*distance))
  Ang = Rad*(180/pi)
  return(Ang)
}

optodat <- read.csv(text="RT,MT,PS,TPS,PAX,TPAX,mx,my,COUNT")
iddat <- read.csv(text = 'PPT,SIDE,VIEW,TRIAL,POSITION,EYE_MOVE,FILENUM')

# reading in all idfiles and compiling
for(i in idfiles){
    tmp <- read.csv(i, sep='\t')[, c(1:6,9)]
    iddat <- rbind(iddat, tmp)
}

row_x <-1

#reading in all TRJ files
for(file in files) {
  tmp <- as.numeric(t(as.vector(read_tsv(file, col_names = FALSE)[1:14, 7])))
  
  optodat[row_x, 1:8] <- tmp[c(1:6,9:10)]
  optodat$COUNT[row_x] <- row_x
  row_x <- row_x + 1
  
}

#finding and removing particular rows from iddat (trials removed from kinematic)
iddat <- iddat[-c(1, 3485, 3895, 3917, 4160),]

#merging id data and optotrak data
res <- merge(data.frame(optodat, row.names = NULL),
             data.frame(iddat, row.names = NULL), by = 0, all = TRUE)[-1] 
res <- res[order(res$COUNT, decreasing = FALSE), ] #ordering by count (so all in order)
# removing DF
res <- res[!res$PPT == 'DMT300', ]
# saving combined file
setwd(anaPath)
write.csv(res, "radial-reaching_allData.csv", row.names = FALSE)

##### calibration data #####
# calibration data-frame
caldat <- res[res$VIEW == 'CAL', ]
colnames(caldat)[colnames(caldat) == 'mx'] <- 'calx' #renaming for merging
colnames(caldat)[colnames(caldat) == 'my'] <- 'caly'
#remove unnecessary trials
caldat <- caldat[c(7,8,10,14)]

res <- res[!res$VIEW =='CAL', ] #data-frame without cal trials, not needed

# merging data-frames with new cal-value
res <- merge(res, caldat, by = c('PPT','POSITION'), all = TRUE) #SO HAPPY THIS FUNCT EXISTS
res[res == -32768] <- NA

##### data analysis #####
# subtracting cal from reach endpoint
res$LANDx <- res$mx - res$calx
res$LANDy <- res$my - res$caly

# counting eye-move and removing eye-move + void
nEye_move <- aggregate(res$EYE_MOVE == '1', by=list(subject_nr = res$PPT), FUN=sum)
nVoid <- aggregate(res$EYE_MOVE == '1', by=list(subject_nr = res$PPT), FUN=sum)
# removing
res <- res[res$EYE_MOVE == 0, c(1:6,9:10,12:13,17:20)]


# plotting to get a look at data

ggplot(res) + geom_point(aes(x = calx, y = caly, colour = POSITION), shape = 3) +
  geom_point(aes(x = mx, y = my, colour = POSITION)) +
  facet_wrap(~PPT*VIEW)
ggsave('radial-reach_Err.png', plot = last_plot(), device = NULL, 
       path = anaPath, scale = 1, width = 15, height = 10, units = 'in')


# tranforming to degrees
res$LANDx_deg <- visAngle(size= res$LANDx, distance = 500) #using visual angle function above
res$LANDy_deg <- visAngle(size= res$LANDy, distance = 500)
# absolute error
res$AE <- sqrt(res$LANDx^2 + res$LANDy^2) #mm
res$AEdeg <- sqrt(res$LANDx_deg^2 + res$LANDy_deg^2) #deg
res$GRP <- factor(substr(res$PPT, 4, 4))

#save compiled data-set
write.csv(res, "radial-reaching_compiled.csv", row.names = FALSE)

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
