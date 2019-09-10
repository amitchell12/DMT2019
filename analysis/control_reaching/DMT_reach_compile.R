library(readr)

mypath <- getwd()

files_201 <- list.files(path=mypath, pattern = "*.TRJ",
                        full.names = TRUE, recursive = TRUE)

optodat <- read.csv(text="SIDE,VIEW,TRIAL,RT,MT,PS,TPS,mx,my")
caldat <- read.csv(text="cx,cy")

row_x <- 1
row_c <- 1

for(i in files_201){
  tmp <- as.numeric(t(as.vector(read_tsv(i, col_names = FALSE)[1:14, 7])))
  SIDE <- substr(i,nchar(i)- 10,nchar(i)- 10)
  VIEW <- substr(i,nchar(i)- 8,nchar(i)- 8)
  TRIAL <- as.numeric(substr(i,nchar(i)- 6,nchar(i)- 4))
  if(VIEW=="F"){VIEW <- "FREE"}else{VIEW<-"PERIPH"}
  if(SIDE == "L"){
  optodat[row_x, 1:9] <- c("LEFT", VIEW, TRIAL, tmp[c(1:4,9:10)])
  row_x <- row_x+1
  } else if (SIDE == "R") {
    optodat[row_x, 1:9] <- c("RIGHT", VIEW, TRIAL, tmp[c(1:4,9:10)])
    row_x <- row_x+1
  } else {
    caldat[row_c, 1:2] <- tmp[10:11]
    row_c <- row_c+1
  }
}

caldat$POSITION <- c(-100,-200,-300,-400,100,200,300,400)

idfiles <- list.files(path=mypath, pattern = "*.dmt",
          full.names = TRUE, recursive = TRUE)

iddat <- read.csv(text="PPT,SIDE,VIEW,TRIAL,POSITION,EYE_MOVE,DATE,TIME")

for(i in idfiles){
  tmp <- read.csv(i)
  iddat <- rbind(iddat, tmp)
}

optodat$TRIAL <- as.numeric(optodat$TRIAL)
optodat$SIDE <- factor(optodat$SIDE)
optodat$VIEW <- factor(optodat$VIEW)

tst <- merge(optodat, iddat)

tst <- merge(tst,caldat)

DMT201 <- tst

DMT201[DMT201 == -32768] <- NA

#DMT201[6] <- as.numeric(DMT201[6])





write.csv(DMT201, "DMT201.csv", row.names = FALSE)

