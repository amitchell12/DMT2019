library(readr)

#on mac
anaPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/analysis/radial_reaching'
dataPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/data'
#on pc
#dataPath <- 'S:/groups/DMT/data'
#anaPath <- 'S:/groups/DMT/analysis/radial_reaching'
setwd(dataPath)

files <- list.files(path=dataPath, pattern = "*.TRJ", full.names = TRUE, recursive = TRUE)
idfiles <- list.files(path=dataPath, pattern = "*.dmt", full.names = TRUE, recursive = TRUE)


optodat <- read.csv(text="RT,MT,PS,TPS,PAX,TPAX,mx,my")
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
  row_x <- row_x + 1
}
  
#finding and removing particular rows from iddat (trials removed from kinematic)


tst <- merge(optodat, iddat)

tst <- merge(tst,caldat)

DMT201 <- tst

DMT201[DMT201 == -32768] <- NA

#DMT201[6] <- as.numeric(DMT201[6])





write.csv(DMT201, "DMT201.csv", row.names = FALSE)

