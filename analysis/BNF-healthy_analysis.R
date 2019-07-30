library(ggplot2)
library(Rmisc)

dataPath ="\\\\chss.datastore.ed.ac.uk/chss/ppls/users/amitch17/Alex_Files/Experiments/DMT2019/DMT2019_rawData/healthy"
setwd(dataPath) #setting the path

library(readr)
#importing correct files
# RIGHT HAND SIDE FIRST
#closed loop beep
h02_clb_right <- 
  read_csv("subject902_pointingTask_CLbeep_right_2019-07-23T15_20_14.csv")
h03_clb_right <-
  read_csv("subject903_pointingTask_CLbeep_right_2019-07-30T10_20_16.csv")
h04_clb_right <-
  read_csv("subject904_pointingTask_CLbeep_right_2019-07-24T12_02_07.csv")
#closed loop norm
h02_cln_right <- 
  read_csv("subject902_pointingTask_CLnorm_right_2019-07-23T16_09_28.csv")
h03_cln_right <-
  read_csv("subject903_pointingTask_CLnorm_right_2019-07-23T16_43_20.csv")
h04_cln_right <-
  read_csv("subject904_pointingTask_CLnorm_right_2019-07-24T11_54_15.csv")
#closed loop dynamic fix
h02_clf_right <- 
  read_csv("subject902_pointingTask_CLfix_right_2019-07-23T16_13_12.csv")
h03_clf_right <-
  read_csv("subject903_pointingTask_CLfix_right_2019-07-23T16_41_44.csv")
h04_clf_right <-
  read_csv("subject904_pointingTask_CLfix_right_2019-07-24T11_58_31.csv")
#open loop beep
h02_ol_right <- 
  read_csv("subject902_pointingTask_OLnorm_right_2019-07-23T16_19_36.csv")
h03_ol_right <-
  read_csv("subject903_pointingTask_OLnorm_right_2019-07-23T16_46_40.csv")
h04_ol_right <-
  read_csv("subject904_pointingTask_OLnorm_right_2019-07-24T12_04_48.csv")

#adding subject number column
h02_clb_right$SUB <- 902
h03_clb_right$SUB <- 903
h04_clb_right$SUB <- 904

h02_cln_right$SUB <- 902
h03_cln_right$SUB <- 903
h04_cln_right$SUB <- 904

h02_clf_right$SUB <- 902
h03_clf_right$SUB <- 903
h04_clf_right$SUB <- 904

h02_ol_right$SUB <- 902
h03_ol_right$SUB <- 903
h04_ol_right$SUB <- 904


#adding type of task
#cl beep
h02_clb_right$TASK <- 1
h03_clb_right$TASK <- 1
h04_clb_right$TASK <- 1
#cl norm
h02_cln_right$TASK <- 2
h03_cln_right$TASK <- 2
h04_cln_right$TASK <- 2
#cl fix
h02_clf_right$TASK <- 3
h03_clf_right$TASK <- 3
h04_clf_right$TASK <- 3
#ol
h02_ol_right$TASK <- 4
h03_ol_right$TASK <- 4
h04_ol_right$TASK <- 4

CLB <- rbind(h02_clb_right, h03_clb_right, h04_clb_right)
CLN <- rbind(h02_cln_right, h03_cln_right, h04_cln_right)
CLF <- rbind(h02_clf_right, h03_clf_right, h04_clf_right)
OL <- rbind(h02_ol_right, h03_ol_right, h04_ol_right)

B <- CLB[c(5,6,17,18,19,20)]
N <- CLN[c(4,5,17,18,19,20)]
F <- CLF[c(7,8,14,15,20,21)]
O <- OL[c(5,6,17,18,19,20)]

names(B)[1:6] <- c("rX", "rY","tX", "tY", "SUB", "TASK")
names(N)[1:6] <- c("rX", "rY","tX", "tY", "SUB", "TASK")
names(F)[1:6] <- c("rX", "rY","tX", "tY", "SUB", "TASK")
names(O)[1:6] <- c("rX", "rY","tX", "tY", "SUB", "TASK")


BNF <- rbind(B, N, F, O)


BNF$COND <- "CLB"
BNF[BNF$TASK == 2, "COND"] <- "CLN"
BNF[BNF$TASK == 3, "COND"] <- "CLF"
BNF[BNF$TASK == 4, "COND"] <- "OL"

BNF$targ <- factor(BNF$tX)

ggplot(BNF)+ geom_point(aes(x=rX, y=rY, colour=targ), size=3) +
              geom_point(aes(x=tX, y = tY, colour=targ), shape=4, size=3, alpha=1) +
              lims(x=c(-50, 900)) + lims(y=c(-300,300)) +
               facet_grid(COND~SUB) +
              theme(axis.text=element_text(size=12),
                       axis.title=element_text(size=12)) +
              theme_bw() -> BNFplot
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
ol_mean = colMeans(O)
ol_error = ol_mean[1]-ol_mean[3]

## NOW THE SAME ON THE LEFT
#closed loop beep
h02_clb_left <- 
  read_csv("subject902_pointingTask_CLbeep_left_2019-07-23T15_31_29.csv")
h03_clb_left <-
  read_csv("subject903_pointingTask_CLbeep_left_2019-07-30T10_21_41.csv")
h04_clb_left <-
  read_csv("subject903_pointingTask_CLfix_left_2019-07-23T16_50_16.csv")
##continue changing here
#closed loop norm
h02_cln_right <- 
  read_csv("subject902_pointingTask_CLnorm_right_2019-07-23T16_09_28.csv")
h03_cln_right <-
  read_csv("subject903_pointingTask_CLnorm_right_2019-07-23T16_43_20.csv")
h04_cln_right <-
  read_csv("subject904_pointingTask_CLnorm_right_2019-07-24T11_54_15.csv")
#closed loop dynamic fix
h02_clf_right <- 
  read_csv("subject902_pointingTask_CLfix_right_2019-07-23T16_13_12.csv")
h03_clf_right <-
  read_csv("subject903_pointingTask_CLfix_right_2019-07-23T16_41_44.csv")
h04_clf_right <-
  read_csv("subject904_pointingTask_CLfix_right_2019-07-24T11_58_31.csv")
#open loop beep
h02_ol_right <- 
  read_csv("subject902_pointingTask_OLnorm_right_2019-07-23T16_19_36.csv")
h03_ol_right <-
  read_csv("subject903_pointingTask_OLnorm_right_2019-07-23T16_46_40.csv")
h04_ol_right <-
  read_csv("subject904_pointingTask_OLnorm_right_2019-07-24T12_04_48.csv")

#adding subject number column
h02_clb_right$SUB <- 902
h03_clb_right$SUB <- 903
h04_clb_right$SUB <- 904

h02_cln_right$SUB <- 902
h03_cln_right$SUB <- 903
h04_cln_right$SUB <- 904

h02_clf_right$SUB <- 902
h03_clf_right$SUB <- 903
h04_clf_right$SUB <- 904

h02_ol_right$SUB <- 902
h03_ol_right$SUB <- 903
h04_ol_right$SUB <- 904


#adding type of task
#cl beep
h02_clb_right$TASK <- 1
h03_clb_right$TASK <- 1
h04_clb_right$TASK <- 1
#cl norm
h02_cln_right$TASK <- 2
h03_cln_right$TASK <- 2
h04_cln_right$TASK <- 2
#cl fix
h02_clf_right$TASK <- 3
h03_clf_right$TASK <- 3
h04_clf_right$TASK <- 3
#ol
h02_ol_right$TASK <- 4
h03_ol_right$TASK <- 4
h04_ol_right$TASK <- 4

CLB <- rbind(h02_clb_right, h03_clb_right, h04_clb_right)
CLN <- rbind(h02_cln_right, h03_cln_right, h04_cln_right)
CLF <- rbind(h02_clf_right, h03_clf_right, h04_clf_right)
OL <- rbind(h02_ol_right, h03_ol_right, h04_ol_right)

B <- CLB[c(5,6,17,18,19,20)]
N <- CLN[c(4,5,17,18,19,20)]
F <- CLF[c(7,8,14,15,20,21)]
O <- OL[c(5,6,17,18,19,20)]

names(B)[1:6] <- c("rX", "rY","tX", "tY", "SUB", "TASK")
names(N)[1:6] <- c("rX", "rY","tX", "tY", "SUB", "TASK")
names(F)[1:6] <- c("rX", "rY","tX", "tY", "SUB", "TASK")
names(O)[1:6] <- c("rX", "rY","tX", "tY", "SUB", "TASK")


BNF <- rbind(B, N, F, O)


BNF$COND <- "CLB"
BNF[BNF$TASK == 2, "COND"] <- "CLN"
BNF[BNF$TASK == 3, "COND"] <- "CLF"
BNF[BNF$TASK == 4, "COND"] <- "OL"

BNF$targ <- factor(BNF$tX)

ggplot(BNF)+ geom_point(aes(x=rX, y=rY, colour=targ), size=3) +
  geom_point(aes(x=tX, y = tY, colour=targ), shape=4, size=3, alpha=1) +
  lims(x=c(-50, 950)) + lims(y=c(-350,350)) +
  facet_grid(COND~SUB) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12)) +
  theme_bw() -> BNFplot
BNFplot
dev.off()