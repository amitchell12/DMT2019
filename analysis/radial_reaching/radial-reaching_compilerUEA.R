library(readr)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(Rmisc)

dataPath <- '/Users/Alex/Documents/DMT/norwich_movement_data'
setwd(dataPath)

files <- list.files(path=dataPath, pattern = "*.TRJ", full.names = TRUE, recursive = TRUE)
idfiles <- list.files(path=dataPath, pattern = "*.csv", full.names = TRUE, recursive = TRUE)

# data-frame structures
optodat <- read.csv(text="RT,MT,PS,TPS,PAX,TPAX,mx,my,FILE")
caldat <- read.csv(text="RT,MT,PS,TPS,PAX,TPAX,mx,my,FILE")
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

    optodat[row_x, 1:9] <- tmp[c(1:6,9:10,15)]
    optodat$COUNT[row_x] <- row_x #counter
    row_x <- row_x + 1
  }
}
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

colnames(caldat)[colnames(caldat) == 'mx'] <- 'calx' #renaming for merging
colnames(caldat)[colnames(caldat) == 'my'] <- 'caly'

# extracting key features from FILE in optodat in the same manner as iddat - for merging
optodat$PPT <- substr(optodat$FILE, 5, 7)
optodat$SIDE <- substr(optodat$FILE, (nchar(optodat$FILE)-15), (nchar(optodat$FILE)-14))
optodat$VIEW <- substr(optodat$FILE, 18, 18)
optodat$TRIAL <- substr(optodat$FILE, (nchar(optodat$FILE)-5), (nchar(optodat$FILE)-4))
# removing '0' that pre-fix digits
optodat$TRIAL <- gsub("(^|[^0-9])0+", "\\1", optodat$TRIAL, perl = TRUE)

# getting 'position' and 'ppt' data for caldat
caldat$PPT <- substr(caldat$FILE, 5, 7)
# position information- from counter and calx
# first remove 6+7 (button presses)
caldat <- subset(caldat, caldat$COUNT !=5 & caldat$COUNT !=6)
# flag position - should just be list 1:4
caldat$POSITION <- c(1:4)
# flag side based on -ve and +ve signed
caldat$SIDE <- factor(caldat$calx < 0, labels = c('Right','Left'))
# extract key variables
caldat <- caldat[, c(11:13,7:8)]

### MERGE
### IDDAT AND OPTODAT TRIALS DO NOT ALIGN - THIS IS AN ISSUE, MAYBE SPEAK TO THE OTHERS ABOUT THIS?
# combining id data and optotrak data
res <- merge(iddat, optodat, by = c('PPT', 'SIDE', 'VIEW', 'TRIAL'))
res <- res[order(res$COUNT, decreasing = FALSE), ] #ordering by count (so all in order)
# renaming conditions to be more sensible :)
res$SIDE <- factor(res$SIDE, labels = c('Left', 'Right'))
res$VIEW <- factor(res$VIEW, labels = c('Free', 'Peripheral'))
# remove file column
res <- res[, c(1:17,19)]


## merging caldat to account for calibrationd data
# labelling 'position' and 'side' in caldat
res <- merge(res, caldat, by = c('PPT','POSITION','SIDE'), all = TRUE) #SO HAPPY THIS FUNCT EXISTS
res <- res[order(res$COUNT, decreasing = FALSE), ] #ordering by count (so all in order)
res[res == -32768] <- NA
res <- res[res$RT > 0 ,] #removing NA values

# counting eye-move and removing eye-move + void
nEye_move <- aggregate(res$RESP == '69', by=list(subject_nr = res$PPT), FUN=sum)
nTimeOut <- aggregate(res$TIMEOUT == '1', by=list(subject_nr = res$PPT), FUN=sum)
nVoid <- aggregate(res$RESP == '73', by=list(subject_nr = res$PPT), FUN=sum)
# removing
res <- res[res$RESP == 32 ,]
res <- res[complete.cases(res) ,] #removing NA values

##### data analysis #####
# subtracting cal from reach endpoint
res$calx <- as.numeric(res$calx)
res$caly <- as.numeric(res$caly)
res$mx <- as.numeric(res$mx)
res$my <- as.numeric(res$my)
res$LANDx <- res$mx - res$calx
res$LANDy <- res$my - res$caly

# plotting to get a look at data - check alignment
res$POSITION <- factor(res$POSITION)
ggplot(res) + geom_point(aes(x = calx, y = caly, colour = POSITION), shape = 3) +
  geom_point(aes(x = mx, y = my, colour = POSITION)) +
  facet_wrap(~PPT+VIEW) +
  ylim(250,500) 
ggsave('radial-reach_Err.png', plot = last_plot(), device = NULL, 
       path = dataPath, scale = 1, width = 15, height = 10, units = 'in')

# absolute error
res$AE <- sqrt(res$LANDx^2 + res$LANDy^2) #mm
res$GRP <- factor(substr(res$PPT, 1, 1))
res$SITE <- 'UEA'

# adding demographic information
# add demographic information to this data
patient_demos <- read.csv('patient_demographics.csv') #loading patient demographics
control_demos <- read.csv('control_demographics.csv') #loading control demos
#extracting ACE data into seperate data-frame
ACEscores <- patient_demos[ ,c(1, 8:13)]
#isolating patient demographic information to bind with control
patient_demos <- patient_demos[, c(1:6)]
demo <- rbind(control_demos, patient_demos)
# changing 'subject_nr' to PPT to match data-set
names(demo)[1] <- 'PPT'

#merging demo with res medians
res <- merge(res, demo, by = 'PPT')

#save compiled data-set
write.csv(res, "radial-reaching_compiled.csv", row.names = FALSE)
