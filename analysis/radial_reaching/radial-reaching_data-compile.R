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
res$POSITION <- factor(res$POSITION)
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
