library(readr)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(Rmisc)
library(dplyr, warn.conflicts = FALSE)
library(forcats)

###### UOE DATA ######
#on mac desktop
#anaPath <- '/Users/Alex/Documents/DMT/analysis/radial_reaching/'
#dataPath <- '/Users/Alex/Documents/DMT/data/'

# work computer
anaPath <- 'S:/groups/DMT/analysis/radial_reaching'
dataPath <- 'S:/groups/DMT/data'
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

optodat <- read.csv(text="px,py,PS,TPS,FILE,COUNT")
iddat <- read.csv(text = 'PPT,SIDE,VIEW,TRIAL,POSITION,EYE_MOVE,FILENUM')
caldat <- read.csv(text='RT,MT,PS,TPS,PAX,TPAX,calx,caly,FILE')

# reading in all idfiles and compiling
for(i in idfiles){
    tmp <- read.csv(i, sep='\t')[, c(1:6,9)]
    iddat <- rbind(iddat, tmp)
}

row_x <-1

#reading in all TRJ files
for(file in files) {
  ## MADE THESE CHANGES - CHECK TOMORROW THAT THEY WORK
  FILE <- basename(file)
  tmp <- c(as.numeric(unlist(read_tsv(file, col_names = FALSE)[1:14, 7])), FILE)
  TPSrow <- (as.numeric(tmp[4])/10)+1
  tmp1 <- as.numeric(t(as.vector(read_tsv(file, col_names = FALSE)[TPSrow, 1:2])))
  
  
  optodat[row_x, 1:2] <- tmp1
  optodat[row_x, 3:5] <- tmp[c(3:4,15)]
  optodat$COUNT[row_x] <- row_x
  row_x <- row_x + 1
  
}

optodat$PPT <- substr(optodat$FILE, 1, 6)
optodat$FILENUM <- substr(optodat$FILE, 8, 10)
optodat$FILENUM <- gsub("(^|[^0-9])0+", "\\1", optodat$FILENUM, perl = TRUE)

## MIGHT NEED TO MAKE CHANGES BELOW - USE UEA IF NEEDED

#merging id data and optotrak data
res <- merge(iddat, optodat, by = c('FILENUM','PPT')) 
res <- res[order(res$COUNT, decreasing = FALSE), ] #ordering by count (so all in order)
# removing DF
res <- res[!res$PPT == 'DMT300', ]
res <- res[!res$VIEW =='CAL', ] #data-frame without cal trials, not needed

##### calibration data #####
# take seperately - different cal data to x,y data here
cal_x <- 1

for(file in files) {
  ## MADE THESE CHANGES - CHECK TOMORROW THAT THEY WORK
  FILE <- basename(file)
  tmp <- c(as.numeric(unlist(read_tsv(file, col_names = FALSE)[1:14, 7])), FILE)
  caldat[cal_x, 1:9] <- tmp[c(1:6,9:10,15)]
  cal_x <- cal_x+1 #counter

}

caldat$PPT <- substr(caldat$FILE, 1, 6)
caldat$FILENUM <- substr(caldat$FILE, 8, 10)
caldat$FILENUM <- gsub("(^|[^0-9])0+", "\\1", caldat$FILENUM, perl = TRUE)

cal <- merge(iddat, caldat, by = c('FILENUM','PPT'))
cal <- cal[cal$VIEW =='CAL', ]
cal <- cal[order(cal$PPT, decreasing = FALSE), ]
# remove DMT300 and extra cal trial on DMT101
cal <- cal[!cal$PPT == 'DMT300', ]
cal <- cal[-c(1) ,]
# remove unnecessary cols
cal <- cal[, c(1:3,6,14,15)]

# merging data-frames with new cal-value
resUOE <- merge(res, cal, by = c('PPT','POSITION','SIDE'), all = TRUE) #SO HAPPY THIS FUNCT EXISTS
resUOE[resUOE == -32768] <- NA
# removing eye-move data
resUOE <- resUOE[resUOE$EYE_MOVE == 0 ,]
resUOE <- resUOE[, c(1:3,5:13,15,16)]


###### UEA DATA ######
# work computer
dataPath <- 'S:/groups/DMT/norwich_movement_data'
setwd(dataPath)

files <- list.files(path=dataPath, pattern = "*.TRJ", full.names = TRUE, recursive = TRUE)
idfiles <- list.files(path=dataPath, pattern = "*.csv", full.names = TRUE, recursive = TRUE)

# data-frame structures
optodat <- read.csv(text="px,py,PS,TPS,FILE,COUNT")
caldat <- read.csv(text="RT,MT,PS,TPS,PAX,TPAX,calx,caly,FILE")
iddat <- read.csv(text = 'POSITION,DELAY,RT,TO,RESP,FILE,TRIAL')

# reading in all idfiles for trial information, and compiling
# excluding practice or voided blocks
for(i in idfiles){
  if (isTRUE(substr(i, nchar(i)-5, nchar(i)-4) == "et")){
    tmp <- read.csv(i, header = FALSE)[, c(1:5)]
    colnames(tmp) <- NULL #removing column names
    tmp$FILE <- basename(i)
    tmp$TRIAL <- 1:length(tmp$FILE)
    iddat <- rbind(iddat, tmp) #data-frame of trial information
  }
}
# adding key trial information
iddat$PPT <- substr(iddat$FILE, 1, 3)
iddat$SIDE <- substr(iddat$FILE, (nchar(iddat$FILE)-11), (nchar(iddat$FILE)-10))
iddat$VIEW <- substr(iddat$FILE, 14, 14)

# re-ordering and naming columns
names(iddat)[1] <- 'POSITION'
names(iddat)[2] <- 'DELAY'
names(iddat)[3] <- 'RTID'
names(iddat)[4] <- 'TIMEOUT'
names(iddat)[5] <- 'RESP'
# reorder so key info in first few columns
iddat <- iddat[, c(8:10,1:5,7)]

row_x <-1
# now for movement (TRJ files), reading in and compiling
for(file in files) {
  # get movement (reaaching data-files)
  if (isTRUE(substr(file, nchar(file)-9, nchar(file)-8) == "et")){
    FILE <- basename(file)
    tmp <- c(as.numeric(unlist(read_tsv(file, col_names = FALSE)[1:14, 7])), FILE)
    TPSrow <- (as.numeric(tmp[4])/5.58) +1
    TPSrow <- round(TPSrow)
    tmp1 <- as.numeric(t(as.vector(read_tsv(file, col_names = FALSE)[TPSrow, 1:2])))
    
    optodat[row_x, 1:2] <- tmp1
    optodat[row_x, 3:5] <- tmp[c(3:4,15)]
    optodat$COUNT[row_x] <- row_x #counter
    row_x <- row_x + 1
  }
}

# extracting key features from FILE in optodat in the same manner as iddat - for merging
optodat$PPT <- substr(optodat$FILE, 5, 7)
optodat$SIDE <- substr(optodat$FILE, (nchar(optodat$FILE)-15), (nchar(optodat$FILE)-14))
optodat$VIEW <- substr(optodat$FILE, 18, 18)
optodat$TRIAL <- substr(optodat$FILE, (nchar(optodat$FILE)-5), (nchar(optodat$FILE)-4))
# removing '0' that pre-fix digits
optodat$TRIAL <- gsub("(^|[^0-9])0+", "\\1", optodat$TRIAL, perl = TRUE)

### calibration data ###
row_x <- 1
cal_x <- 1
# get calibration files
for(file in files) {
  # get movement (reaaching data-files)
  if (isTRUE(substr(file, nchar(file)-13, nchar(file)-12) == "st")){
    FILE <- basename(file)
    tmp <- c(as.numeric(unlist(read_tsv(file, col_names = FALSE)[1:14, 7])), FILE)
    
    caldat[row_x, 1:9] <- tmp[c(1:6,9:10,15)]
    caldat$COUNT[row_x] <- cal_x #counter for each calib trial
    row_x <- row_x + 1
    if(isTRUE(cal_x >= 10)){
      cal_x <- (cal_x+1)-10
    }
    else{
      cal_x <- cal_x + 1
    }
  }
}


# getting 'position' and 'ppt' data for caldat
caldat$PPT <- substr(caldat$FILE, 5, 7)

# saving button locations for calculating amplitude and angular error
button_left <- caldat[caldat$COUNT == 5 ,]
button_right <- caldat[caldat$COUNT == 6 ,]
# get side
button_left$SIDE <- 'Left'
button_right$SIDE <- 'Right'
buttons <- rbind(button_left, button_right)
buttons$calx <- as.numeric(buttons$calx)
buttons$caly <- as.numeric(buttons$caly)
# one erroneous trial, correct
for (i in 1:length(buttons$PPT)){
  if (isTRUE(buttons$calx[i] < -100)){
    buttons$calx[i] = -30
    buttons$caly[i] = 39
  }
}
# rename columns
names(buttons)[7] <- 'sX'
names(buttons)[8] <- 'sY'
# remove unnecessary data
buttons <- buttons[, c(7,8,11,12)]

# position information- from counter and calx
# first remove 6+7 (button presses)
caldat <- subset(caldat, caldat$COUNT !=5 & caldat$COUNT !=6)
# flag position - should just be list 1:4
caldat$POSITION <- c(1:4)
# flag side based on -ve and +ve signed
caldat$SIDE <- factor(caldat$calx < 0, labels = c('Right','Left'))
# extract key variables
caldat <- caldat[, c(11:13,7:8)]
## merge button data to caldat
caldat <- merge(caldat, buttons, by = c('SIDE','PPT'))

## missing calibration data - load this
cal_fits <- read.csv(text='TRIAL,rX,rY,FILE')
for(i in idfiles){
  if (isTRUE(substr(basename(i), 11, 14) == 'fits')){
    tmp <- read.csv(i, header = FALSE)[, c(1,2,4)]
    tmp$FILE <- basename(i)
    cal_fits <- rbind(cal_fits, tmp) #data-frame of trial information
  }
}

# re-arrage to fit other calibration data, then bind and order
colnames(cal_fits)[colnames(cal_fits) == 'V2'] <- 'calx' #renaming for merging
colnames(cal_fits)[colnames(cal_fits) == 'V4'] <- 'caly'
cal_fits$PPT <- substr(cal_fits$FILE, 1, 3)
cal_fits$POSITION <- c(1:4)
cal_fits$SIDE <- factor(cal_fits$calx < 0, labels = c('Right','Left'))
## add missing button data to cal_fits
# estimate start position from the mean of the other values - for missing data
mean_startX <- summarySEwithin(data = buttons, measurevar = 'sX', withinvars = 'SIDE')
mean_startY <- summarySEwithin(data = buttons, measurevar = 'sY', withinvars = 'SIDE')
# adding mean fitted data to missing data-set
for (i in 1:length(cal_fits$PPT)){
  if (isTRUE(cal_fits$SIDE[i] == 'Left')){
    cal_fits$sX[i] <- mean_startX$sX[1]
    cal_fits$sY[i] <- mean_startY$sY[1]
  }
  else if (isTRUE(cal_fits$SIDE[i] == 'Right')){
    cal_fits$sX[i] <- mean_startX$sX[2]
    cal_fits$sY[i] <- mean_startY$sY[2]
  }
}

# merge 
cal_fits <- cal_fits[, c(5:7,2,3,8,9)]
caldat <- rbind(caldat,cal_fits)

### MERGE
# combining id data and optotrak data
resUEA <- merge(iddat, optodat, by = c('PPT', 'SIDE', 'VIEW', 'TRIAL'))
resUEA <- resUEA[order(resUEA$COUNT, decreasing = FALSE), ] #ordering by count (so all in order)
# renaming conditions to be more sensible :)
resUEA$SIDE <- factor(resUEA$SIDE, labels = c('Left', 'Right'))
resUEA$VIEW <- factor(resUEA$VIEW, labels = c('Free', 'Peripheral'))

# remove file column
resUEA <- resUEA[, c(1:13,15)]

## merging caldat to account for calibrationd data
# labelling 'position' and 'side' in caldat
resUEA <- merge(resUEA, caldat, by = c('PPT','POSITION','SIDE'), all = TRUE) #SO HAPPY THIS FUNCT EXISTS
resUEA <- resUEA[order(resUEA$COUNT, decreasing = FALSE), ] #ordering by count (so all in order)
resUEA[resUEA == -32768] <- NA
# counting eye-move and removing eye-move + void
nEye_move <- aggregate(resUEA$RESP == '69', by=list(subject_nr = resUEA$PPT), FUN=sum)
nTimeOut <- aggregate(resUEA$TIMEOUT == '1', by=list(subject_nr = resUEA$PPT), FUN=sum)
nVoid <- aggregate(resUEA$RESP == '73', by=list(subject_nr = resUEA$PPT), FUN=sum)
# removing
resUEA <- resUEA[resUEA$RESP == 32 ,]
resUEA <- resUEA[complete.cases(resUEA) ,] #removing NA values

##### COMBINING UEA AND UOE DATA ######
# UOE res - to match UEA
resUOE$PPT <- substr(resUOE$PPT, 4, 6)
resUOE$VIEW <- factor(resUOE$VIEW, labels = c('Free', 'Peripheral'))
resUOE$SIDE <- factor(resUOE$SIDE, labels = c('Left', 'Right'))
# adding start positions of finger
resUOE$sX <- 500
resUOE$sY <- 0

# UEA res - match data
# remove participants
resUEA <- resUEA[resUEA$PPT != '311' ,] #participant 311 had TIA
# first 10 participants did different set-up to UEA patients, also remove
oldPP <- c(301:310)
resUEA <- subset(resUEA, ! PPT %in% oldPP)
# change position labelling here too - needs to match UOE (100,200,300,400mm)
# split into two data-frames (left/right sides) and label, the bind again
UEAright <- resUEA[resUEA$SIDE == 'Right' ,]
UEAleft <- resUEA[resUEA$SIDE == 'Left' ,]
# for right 1= closest, 4= furthest away
UEAright$POSITION <- factor(UEAright$POSITION, labels = c('100','200','300','400'))
# for left 1=furthest, 4= nearest
UEAleft$POSITION <- factor(UEAleft$POSITION, labels = c('-400','-300','-200','-100'))
# bind back together!
resUEA <- rbind(UEAright, UEAleft)

# site labelling
resUOE$SITE <- 'UOE'
resUEA$SITE <- 'UEA'

# removing unnecessary cols
resUOE <- resUOE[, c(1:5,7:10,12:17)]
resUEA <- resUEA[, c(1:5,10:19)]

# bind
resALL <- rbind(resUOE, resUEA)
setwd(anaPath)
write.csv(resALL, 'RAD_PSposition.csv', row.names = FALSE)

##### PREP FOR DATA ANALYSIS #####
# getting group and diagnosis
resALL$GRP <- factor(substr(resALL$PPT, 1, 1))
# adding demographic information
setwd('S:/groups/DMT/data')
patient_demos <- read.csv('patient_demographics.csv') #loading patient demographics
control_demos <- read.csv('control_demographics.csv') #loading control demos
#extracting ACE data into seperate data-frame
ACEscores <- patient_demos[ ,c(1, 8:13)]
#isolating patient demographic information to bind with control
patient_demos <- patient_demos[, c(1:6)]
demo <- rbind(control_demos, patient_demos)
# changing 'subject_nr' to PPT to match data-set
names(demo)[1] <- 'PPT'

#merging demo with res
resALL <- merge(demo, resALL, by = 'PPT')
setwd(anaPath)
write.csv(resALL, 'RAD_PSposition.csv', row.names = FALSE)

# subtracting cal from reach endpoint
res$LANDx <- res$px - res$calx
res$LANDy <- res$py - res$caly


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
