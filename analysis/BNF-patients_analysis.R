library(ggplot2)
library(Rmisc)

dataPath ="\\\\chss.datastore.ed.ac.uk/chss/ppls/users/amitch17/Alex_Files/Experiments/DMT2019/DMT2019_rawData/patients"
setwd(dataPath) #setting the path

library(readr)
#importing correct files
#closed loop beep
p01_clb_right <- 
  read_csv("subject801_pointingTask_CLbeep_right_2019-07-26T14_26_13.csv")
p02_clb_right <- 
  read_csv("subject802_pointingTask_CLbeep_right_2019-07-26T14_28_59.csv")
p03_clb_right <-
  read_csv("subject803_pointingTask_CLbeep_right_2019-07-26T14_49_43.csv")
p04_clb_right <-
  read_csv("subject804_pointingTask_CLbeep_right_2019-07-26T14_51_38.csv")
#closed loop norm
p01_cln_right <- 
  read_csv("subject801_pointingTask_CLnorm_right_2019-07-26T14_32_14.csv")
p02_cln_right <- 
  read_csv("subject802_pointingTask_CLnorm_right_2019-07-26T14_32_23.csv")
p03_cln_right <-
  read_csv("subject803_pointingTask_CLnorm_right_2019-07-26T14_57_37.csv")
p04_cln_right <-
  read_csv("subject804_pointingTask_CLnorm_right_2019-07-26T14_56_15.csv")
#closed loop dynamic fix
p01_clf_right <- 
  read_csv("subject801_pointingTask_CLfix_right_2019-07-26T14_37_16.csv")
p02_clf_right <- 
  read_csv("subject802_pointingTask_CLfix_right_2019-07-26T14_37_42.csv")
p03_clf_right <-
  read_csv("subject803_pointingTask_CLfix_right_2019-07-26T15_01_30.csv")
p04_clf_right <-
  read_csv("subject804_pointingTask_CLfix_right_2019-07-26T15_00_04.csv")

#adding subject number column
p01_clb_right$SUB <- 801
p02_clb_right$SUB <- 802
p03_clb_right$SUB <- 803
p04_clb_right$SUB <- 804

p01_cln_right$SUB <- 801
p02_cln_right$SUB <- 802
p03_cln_right$SUB <- 803
p04_cln_right$SUB <- 804

p01_clf_right$SUB <- 801
p02_clf_right$SUB <- 802
p03_clf_right$SUB <- 803
p04_clf_right$SUB <- 804

#adding type of task
p01_clb_right$TASK <- 1 #cl beep
p02_clb_right$TASK <- 1
p03_clb_right$TASK <- 1
p04_clb_right$TASK <- 1

p01_cln_right$TASK <- 2 #cl norm
p02_cln_right$TASK <- 2
p03_cln_right$TASK <- 2
p04_cln_right$TASK <- 2

p01_clf_right$TASK <- 3 #cl fix
p02_clf_right$TASK <- 3
p03_clf_right$TASK <- 3
p04_clf_right$TASK <- 3


CLB <- rbind(p01_clb_right, p02_clb_right, p03_clb_right, p04_clb_right)

CLN <- rbind(p01_cln_right, p02_cln_right, p03_cln_right, p04_cln_right)

CLF <-  rbind(p01_clf_right, p02_clf_right, p03_clf_right, p04_clf_right)

B <- CLB[c(5,6,14,15,19,20)]
N <- CLN[c(4,5,13,14,19,20)]
F <- CLF[c(7,8,14,15,20,21)]

names(B)[1:6] <- c("rX", "rY","tX", "tY", "SUB", "TASK")
names(N)[1:6] <- c("rX", "rY","tX", "tY", "SUB", "TASK")
names(F)[1:6] <- c("rX", "rY","tX", "tY", "SUB", "TASK")


BNF <- rbind(B, N, F)


BNF$PAT <- "MCI"
BNF[BNF$SUB == 801, "PAT"] <- "AD-1"
BNF[BNF$SUB == 802, "PAT"] <- "AD-2"
BNF[BNF$SUB == 803, "PAT"] <- "MCI-1"
BNF[BNF$SUB == 804, "PAT"] <- "MCI-2"

BNF$COND <- "CLB"
BNF[BNF$TASK == 2, "COND"] <- "CLN"
BNF[BNF$TASK == 3, "COND"] <- "CLF"

BNF$targ <- factor(BNF$tX)

ggplot(BNF)+ geom_point(aes(x=rX, y=rY, colour=targ), size=3) +
              geom_point(aes(x=tX, y = tY, colour=targ), shape=4, size=3, alpha=1) +
              lims(x=c(-50, 900)) + lims(y=c(-300,300)) +
               facet_grid(COND~PAT) +
              theme(axis.text=element_text(size=12),
                       axis.title=element_text(size=12)) +
              theme_bw() -> BNFplot

png("BNF-patient.png", width = 30, height =20)
pdf("BNF-patient.pdf",width=30,height=20,paper='special') 
BNFplot
dev.off()

#mean error
library(matrixStats)
clb_mean = colMeans(B)
clb_error = clb_mean[1]-clb_mean[3]
cln_mean = colMeans(N)
cln_error = cln_mean[1]-cln_mean[3]
clf_mean = colMeans(F)
clf_error = clf_mean[1]-clf_mean[3]

