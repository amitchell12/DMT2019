library(readr)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(Rmisc)
library(dplyr, warn.conflicts = FALSE)
library(forcats)

###### UOE DATA ######
#on mac
#anaPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/analysis/radial_reaching'
#dataPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/data'
#on mac desktop
anaPath <- '/Users/Alex/Documents/DMT/analysis/radial_reaching/'
dataPath <- '/Users/Alex/Documents/DMT/data/'
setwd(dataPath)

files <- list.files(path=dataPath, pattern = "*.TRJ", full.names = TRUE, recursive = TRUE)
idfiles <- list.files(path=dataPath, pattern = "*.dmt", full.names = TRUE, recursive = TRUE)

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
  ## MADE THESE CHANGES - CHECK TOMORROW THAT THEY WORK
  tmp <- as.numeric(t(as.vector(read_tsv(file, col_names = FALSE)[1:14, 7])))
  TPSrow <- (tmp[4]/10)+1
  tmp1 <- as.numeric(t(as.vector(read_tsv(file, col_names = FALSE)[TPSrow, 1:3])))
  
  optodat[row_x, 1:3] <- tmp1
  optodat$COUNT[row_x] <- row_x
  row_x <- row_x + 1
  
}

## MIGHT NEED TO MAKE CHANGES BELOW - USE UEA IF NEEDED

#merging id data and optotrak data
res <- merge(data.frame(optodat, row.names = NULL),
             data.frame(iddat, row.names = NULL), by = 0, all = TRUE)[-1] 
res <- res[order(res$COUNT, decreasing = FALSE), ] #ordering by count (so all in order)
# removing DF
res <- res[!res$PPT == 'DMT300', ]

##### calibration data #####
# calibration data-frame
caldat <- res[res$VIEW == 'CAL', ]
colnames(caldat)[colnames(caldat) == 'mx'] <- 'calx' #renaming for merging
colnames(caldat)[colnames(caldat) == 'my'] <- 'caly'
#remove unnecessary trials
caldat <- caldat[c(7,8,10,14)]
#remove first row of caldat data - incorrect calib
caldat <- caldat[-c(1) ,]

res <- res[!res$VIEW =='CAL', ] #data-frame without cal trials, not needed

# merging data-frames with new cal-value
res <- merge(res, caldat, by = c('PPT','POSITION'), all = TRUE) #SO HAPPY THIS FUNCT EXISTS
res[res == -32768] <- NA

##### PREP FOR DATA ANALYSIS #####
# subtracting cal from reach endpoint
res$LANDx <- res$mx - res$calx
res$LANDy <- res$my - res$caly

#renaming some stuff
res$GRP <- factor(substr(res$PPT, 4, 4))
res$PPT <- substr(res$PPT, 4, 6)
res$SITE <- 'UOE'

# plotting to get a look at data
res$POSITION <- factor(res$POSITION)
ggplot(res) + geom_point(aes(x = calx, y = caly, colour = POSITION), shape = 3) +
  geom_point(aes(x = mx, y = my, colour = POSITION)) +
  facet_wrap(~PPT*VIEW)
ggsave('radial-reach_Err.png', plot = last_plot(), device = NULL, 
       path = anaPath, scale = 1, width = 15, height = 10, units = 'in')

## data transformations and calculations
# tranforming to degrees
res$LANDx_deg <- visAngle(size= res$LANDx, distance = 500) #using visual angle function above
res$LANDy_deg <- visAngle(size= res$LANDy, distance = 500)

# absolute error - in mm
res$AE <- sqrt(res$LANDx^2 + res$LANDy^2) #mm

## start position
res$sX <- 500
res$sY <- 0

## angular error and amplitude error
# calculating target and end-point position relative to start-point
res$rX <- res$mx - res$sX
res$rY <- res$my - res$sY
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

# adding demographic information
patient_demos <- read.csv('patient_demographics.csv') #loading patient demographics
control_demos <- read.csv('control_demographics 2.csv') #loading control demos
#extracting ACE data into seperate data-frame
ACEscores <- patient_demos[ ,c(1, 8:13)]
#isolating patient demographic information to bind with control
patient_demos <- patient_demos[, c(1:6)]
demo <- rbind(control_demos, patient_demos)
# changing 'subject_nr' to PPT to match data-set
names(demo)[1] <- 'PPT'

#merging demo with res
res <- merge(demo, res, by = 'PPT')

# counting eye-move and removing eye-move + void
nEye_move <- aggregate(EYE_MOVE == 1 ~ PPT * diagnosis * VIEW, sum, data = res)
nVoid <- aggregate(EYE_MOVE == -1 ~ diagnosis * VIEW, sum, data = res)

nTrials <- res %>% 
  group_by(PPT,VIEW) %>% 
  tally()
nTrials <- na.omit(nTrials)

# getting percentage eye-move for each group, and each condition
Eye_move <- merge(nEye_move, nTrials)
names(Eye_move)[4] <- 'eye_move'
sumEye_move <- aggregate(eye_move ~ diagnosis * VIEW, sum, data = Eye_move)
sumTrials <- aggregate(n ~ diagnosis * VIEW, sum, data = Eye_move)
totEye_move <- merge(sumEye_move, sumTrials)

# save this to combine with UEA data
setwd(anaPath)
write.csv(totEye_move, 'edinburgh_eye-move.csv', row.names = FALSE)

# removing
res <- res[res$EYE_MOVE == 0, c(1:7,17:18,8:15,22:42)]

#do some pre-emptive renaming
res$VIEW <- factor(res$VIEW, labels = c('Free', 'Peripheral'))
res$SIDE <- factor(res$SIDE, labels = c('Left','Right'))

#save compiled data-set
write.csv(res, "radial-reaching_compiled.csv", row.names = FALSE)
