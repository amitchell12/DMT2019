library(ggplot2)
library(Rmisc)

dataPath = "/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/DMT2019_rawdata/patients/"
setwd(dataPath) #setting the path

library(readr)
subject801_pointingTask_CLbeep_right <- 
  read_csv("subject801_pointingTask_CLbeep_right_2019-07-26T14_26_13.csv")

CLB <- rbind(subject801_pointingTask_CLbeep_right_2019.07.26T14_26_13,
      subject802_pointingTask_CLbeep_right_2019.07.26T14_49_43)

CLN <- rbind(subject801_pointingTask_CLnorm_right_2019.07.26T14_32_14,
             subject802_pointingTask_CLnorm_right_2019.07.26T14_57_37)

CLF <- rbind(subject801_pointingTask_CLfix_right_2019.07.26T14_37_16,
             subject802_pointingTask_CLfix_right_2019.07.26T15_01_30)

B <- CLB[c(5,6,14,15,19)]
N <- CLN[c(4,5,13,14,19)]
F <- CLF[c(7,8,14,15,20)]

names(B)[1:5] <- c("rX", "rY","tX", "tY", "SUB")
names(N)[1:5] <- c("rX", "rY","tX", "tY", "SUB")
names(F)[1:5] <- c("rX", "rY","tX", "tY", "SUB")


BNF <- rbind(B, N, F)


BNF$PAT <- "MCI"
BNF[BNF$SUB > 801-3, "PAT"] <- "AD"

BNF$COND <- "CLB"
BNF[BNF$SUB %in% c(801-2,803-2), "COND"] <- "CLN"
BNF[BNF$SUB %in% c(801-3,802-3), "COND"] <- "CLF"

BNF$targ <- factor(BNF$tX)

ggplot(BNF)+ geom_point(aes(x=rX, y=rY, colour=targ), size=3) +
              geom_point(aes(x=tX, y = tY, colour=targ), shape=4, size=3, alpha=1) +
              lims(x=c(0, 800)) + 
              facet_grid(COND~PAT) +
              theme_bw() -> BNFplot

png("BNF.png", res=300, units="in", width = 12, height =9)
BNFplot
dev.off()


