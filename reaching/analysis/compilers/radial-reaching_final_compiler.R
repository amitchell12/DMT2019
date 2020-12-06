## RADIAL REACHING FINAL COMPILER ##
# compiling UEA and UOE data together into one main doc #

# libraries
library(plyr)

#PATHS
anaPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/analysis/radial_reaching'
UEAPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/norwich_movement_data'
dataPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/data'

setwd(anaPath) #for Edinburgh data 
resUOE <- read.csv('radial-reaching_compiledUOE.csv')
setwd(UEAPath) #for Norwich data
resUEA <- read.csv('radial-reaching_compiledUEA.csv') 

##### DATA ORGANISE #####
# UOE data, getting rid of 'deg' values - just use mm
resUOE <- resUOE[, c(1:23,26:38)]
# UEA data getting rid of not needed values
resUEA <- resUEA[, c(1:4,10:17,19:42)]

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

# now let's bind! & order by participant
res <- rbind(resUOE, resUEA)
res <- res[order(res$PPT),]

# changing levelsfor plotting
res$VIEW <- factor(res$VIEW) #changing so only 2 levels recorded
res$GRP <- factor(res$GRP)
#levels(res$GRP) <- c('Control', 'Patient') #changing group name from 1 = control, 2 = AD
res$diagnosis <- factor(res$diagnosis)
res$POSITION <- factor(res$POSITION)
# capitalising names of demographics to match rest
colnames(res)[colnames(res) == 'gender'] <- 'GENDER'
colnames(res)[colnames(res) == 'age'] <- 'AGE'
colnames(res)[colnames(res) == 'hand'] <- 'HAND'
colnames(res)[colnames(res) == 'education'] <- 'ED'
colnames(res)[colnames(res) == 'diagnosis'] <- 'DIAGNOSIS'

## getting dominant + non-dominant sides, for analysis
res$HAND <- factor(res$HAND, labels = c('Left','Right'))
# if hand = side, dominant; else non dominant
res$DOM <- as.numeric(res$SIDE == res$HAND) #1 = dominant, 0 = non-dominant
res$DOM <- factor(res$DOM, labels= c('ND','D'))
# change order so important var up-front
res <- res[, c(1:7,37,16,17,8:15,18:37)]

##### TRIAL OUTLIERS #####
# calculating z-score for each participant
# get PP matrix
PPT <- count(res, vars = PPT)
Z <- read.csv(text = 'AEZ,AE,PPT,VIEW')

# first for free reaching
resFree <- res[res$VIEW == 'Free' ,]
for (p in PPT$vars){
  tmp <- resFree[resFree$PPT == p ,]
  # temporary matrix to save data to
  Zscore <- data.frame(matrix(ncol = 0, nrow = length(tmp$PPT)))
  Zscore$AEZ <- abs(scale(tmp$AE))
  Zscore$AE <- tmp$AE
  Zscore$PPT <- p
  Zscore$VIEW <- 'Free'
  
  Z <- rbind(Z, Zscore)
}
# then for peripheral reaching
resPeriph <- res[res$VIEW == 'Peripheral' ,]
for (p in PPT$vars){
  tmp <- resPeriph[resPeriph$PPT == p ,]
  # temporary matrix to save data to
  Zscore <- data.frame(matrix(ncol = 0, nrow = length(tmp$PPT)))
  Zscore$AEZ <- abs(scale(tmp$AE))
  Zscore$AE <- tmp$AE
  Zscore$PPT <- p
  Zscore$VIEW <- 'Peripheral'
  
  Z <- rbind(Z, Zscore)
}

# merge data-frames
res <- merge(res,Z)
## remove trials where absolute z-score is > 4
XCLUDE <- res[res$AEZ > 4 ,]
res <- res[!(res$AEZ %in% XCLUDE$AEZ), ]

# save Edinburgh & UEA compiled data to ana path
setwd(anaPath)
write.csv(res, 'all_radial-reaching_compiled.csv', row.names = FALSE)

