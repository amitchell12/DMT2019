library(readr)
library(ggplot2)

#on mac
#anaPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/analysis/radial_reaching'
#dataPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/data'
#on pc
dataPath <- 'S:/groups/DMT/data'
anaPath <- 'S:/groups/DMT/analysis/radial_reaching'
setwd(dataPath)

files <- list.files(path=dataPath, pattern = "*.TRJ", full.names = TRUE, recursive = TRUE)
idfiles <- list.files(path=dataPath, pattern = "*.dmt", full.names = TRUE, recursive = TRUE)


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
iddat <- iddat[-c(122, 3485, 3895, 3917, 4160),]

#merging id data and optotrak data
res <- merge(data.frame(optodat, row.names = NULL),
             data.frame(iddat, row.names = NULL), by = 0, all = TRUE)[-1] 
res <- res[order(res$COUNT, decreasing = FALSE), ] #ordering by count (so all in order)
# removing DF
res <- res[!res$PPT == 'DMT300', ]
# saving combined file
setwd(anaPath)
write.csv(res, "radial-reaching_compiled.csv", row.names = FALSE)

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

# plotting to get a look at data
ggplot(res) + geom_point(aes(x = calx, y = caly, colour = POSITION), shape = 3) +
  geom_point(aes(x = mx, y = my, colour = POSITION)) +
  facet_wrap(~PPT*VIEW)

##### data analysis #####
# subtracting cal from reach endpoint
res$LANDx <- res$mx - res$calx
res$LANDy <- res$my - res$caly

# plotting position vs land_y
ggplot(res) + geom_point(aes(x = POSITION, y = LANDy, colour = VIEW)) +
  facet_wrap(~PPT) + ylim(-100, 100)

ggsave('radial-reach_yErr.png', plot = last_plot(), device = NULL, 
       path = anaPath, scale = 1, width = 15, height = 10, units = 'in')




