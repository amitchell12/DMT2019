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
h02_clb_right$TASK <- "clb"
h03_clb_right$TASK <- "clb"
h04_clb_right$TASK <- "clb"
#cl norm
h02_cln_right$TASK <- "cln"
h03_cln_right$TASK <- "cln"
h04_cln_right$TASK <- "cln"
#cl fix
h02_clf_right$TASK <- "clf"
h03_clf_right$TASK <- "clf"
h04_clf_right$TASK <- "clf"
#ol
h02_ol_right$TASK <- "ol"
h03_ol_right$TASK <- "ol"
h04_ol_right$TASK <- "ol"

CLB <- rbind(h02_clb_right, h03_clb_right, h04_clb_right)
CLN <- rbind(h02_cln_right, h03_cln_right, h04_cln_right)
CLF <- rbind(h02_clf_right, h03_clf_right, h04_clf_right)
OL <- rbind(h02_ol_right, h03_ol_right, h04_ol_right)

CLB_r <- CLB[c(5,6,17,18,19,20)]
CLN_r<- CLN[c(4,5,17,18,19,20)]
CLF_r <- CLF[c(7,8,14,15,20,21)]
OL_r <- OL[c(5,6,17,18,19,20)]

names(CLB_r)[1:6] <- c("rX", "rY","tX", "tY", "SUB", "TASK")
names(CLN_r)[1:6] <- c("rX", "rY","tX", "tY", "SUB", "TASK")
names(CLF_r)[1:6] <- c("rX", "rY","tX", "tY", "SUB", "TASK")
names(OL_r)[1:6] <- c("rX", "rY","tX", "tY", "SUB", "TASK")


rData <- rbind(CLB_r,CLN_r,CLF_r,OL_r)

rData$COND <- "CLB"
rData[rData$TASK == "cln", "COND"] <- "CLN"
rData[rData$TASK == "clf", "COND"] <- "CLF"
rData[rData$TASK == "ol", "COND"] <- "OL"

rData$targ <- factor(rData$tX)

ggplot(rData)+ geom_point(aes(x=rX, y=rY, colour=targ), size=3) +
  geom_point(aes(x=tX, y = tY, colour=targ), shape=4, size=3, alpha=1) +
  lims(x=c(-50, 950)) + lims(y=c(-350,350)) +
  facet_grid(COND~SUB) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12)) +
  theme_bw() -> RDplot
RDplot
dev.off()

#mean error
#change this for later 
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
  read_csv("subject904_pointingTask_CLbeep_left_2019-07-24T12_10_50.csv")
#closed loop norm
#continue from here
h02_cln_left <- 
  read_csv("subject902_pointingTask_CLnorm_left_2019-07-23T16_21_49.csv")
h03_cln_left <-
  read_csv("subject903_pointingTask_CLnorm_left_2019-07-23T16_51_52.csv")
h04_cln_left <-
  read_csv("subject904_pointingTask_CLnorm_left_2019-07-24T12_07_22.csv")
#closed loop dynamic fix
h02_clf_left <- 
  read_csv("subject902_pointingTask_CLfix_left_2019-07-23T16_24_06.csv")
h03_clf_left <-
  read_csv("subject903_pointingTask_CLfix_left_2019-07-23T16_50_16.csv")
h04_clf_left <-
  read_csv("subject904_pointingTask_CLfix_left_2019-07-24T12_09_04.csv")
#open loop beep
h02_ol_left <- 
  read_csv("subject902_pointingTask_OLnorm_left_2019-07-23T16_28_17.csv")
h03_ol_left <-
  read_csv("subject903_pointingTask_OLnorm_left_2019-07-23T16_54_59.csv")
h04_ol_left <-
  read_csv("subject904_pointingTask_OLnorm_left_2019-07-24T12_12_30.csv")

#dropping sanity check 
#adding subject number column
h02_clb_left$SUB <- 902
h03_clb_left$SUB <- 903
h04_clb_left$SUB <- 904

h02_cln_left$SUB <- 902
h03_cln_left$SUB <- 903
h04_cln_left$SUB <- 904

h02_clf_left$SUB <- 902
h03_clf_left$SUB <- 903
h04_clf_left$SUB <- 904

h02_ol_left$SUB <- 902
h03_ol_left$SUB <- 903
h04_ol_left$SUB <- 904


#adding type of task
#cl beep
h02_clb_left$TASK <- "clb"
h03_clb_left$TASK <- "clb"
h04_clb_left$TASK <- "clb"
#cl norm
h02_cln_left$TASK <- "cln"
h03_cln_left$TASK <- "cln"
h04_cln_left$TASK <- "cln"
#cl fix
h02_clf_left$TASK <- "clf"
h03_clf_left$TASK <- "clf"
h04_clf_left$TASK <- "clf"
#ol
h02_ol_left$TASK <- "ol"
h03_ol_left$TASK <- "ol"
h04_ol_left$TASK <- "ol"

CLB <- rbind(h02_clb_left, h03_clb_left, h04_clb_left)
CLN <- rbind(h02_cln_left, h03_cln_left, h04_cln_left)
CLF <- rbind(h02_clf_left, h03_clf_left, h04_clf_left)
OL <- rbind(h02_ol_left, h03_ol_left, h04_ol_left)

CLB_l <- CLB[c(5,6,17,18,19,20)]
CLN_l<- CLN[c(4,5,17,18,19,20)]
CLF_l <- CLF[c(7,8,14,15,20,21)]
OL_l <- OL[c(5,6,17,18,19,20)]

names(CLB_l)[1:6] <- c("rX", "rY","tX", "tY", "SUB", "TASK")
names(CLN_l)[1:6] <- c("rX", "rY","tX", "tY", "SUB", "TASK")
names(CLF_l)[1:6] <- c("rX", "rY","tX", "tY", "SUB", "TASK")
names(OL_l)[1:6] <- c("rX", "rY","tX", "tY", "SUB", "TASK")


lData <- rbind(CLB_l, CLN_l, CLF_l, OL_l)

lData$COND <- "CLB"
lData[lData$TASK == "cln", "COND"] <- "CLN"
lData[lData$TASK == "clf", "COND"] <- "CLF"
lData[lData$TASK == "ol", "COND"] <- "OL"

lData$targ <- factor(lData$tX)

ggplot(lData)+ geom_point(aes(x=rX, y=rY, colour=targ), size=3) +
  geom_point(aes(x=tX, y = tY, colour=targ), shape=4, size=3, alpha=1) +
  lims(x=c(-950, 50)) + lims(y=c(-350,350)) +
  facet_grid(COND~SUB) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12)) +
  theme_bw() -> LDplot
LDplot
dev.off()

## EXTRACT SUMMARY DATA ##
#directional error
rData$x_err = rData$rX - rData$tX
lData$x_err = lData$rX - lData$tX
rData$y_err = rData$rY - rData$tY
lData$y_err = lData$rY - lData$tY

# absolute error
rData$abs_err = sqrt(rData$x_err^2 + rData$y_err^2)
lData$abs_err = sqrt(lData$x_err^2 + lData$y_err^2)

#ecccentricity coding
rData$ECC <- factor(cut(rData$tX, 3, label=c("N", "M", "F")))
lData$ECC <- factor(cut(lData$tX, 3, label=c("F", "M", "N"))) #has to be opposite, sign reversed
#making factors
rData$SUB <- factor(rData$SUB)
lData$SUB <- factor(lData$SUB)
rData$TASK <- factor(rData$TASK)
lData$TASK <- factor(lData$TASK)

data <- rbind(rData,lData)
data$side <- "right"
data[data$tX <0, "side"] <- "left"

#plotting
ggplot(rData, aes(x=SUB, y=abs_err, colour = TASK)) +
  geom_boxplot(outlier.alpha=0) +
  geom_jitter(position=position_dodge(width=.8), size=4, alpha=.5) +
  ylim(0,300) + geom_hline(yintercept=55, linetype="dotted") +
  theme_bw() + theme(legend.position = "bottom") -> rDataplot 
rDataplot
#same for left
ggplot(lData, aes(x=SUB, y=abs_err, colour = TASK)) +
  geom_boxplot(outlier.alpha=0) +
  geom_jitter(position=position_dodge(width=.8), size=4, alpha=.5) +
  ylim(0,300) + geom_hline(yintercept=55, linetype="dotted") +
  theme_bw() + theme(legend.position = "bottom") -> lDataplot 
lDataplot
  
library(ggpubr)
taskFig <- ggarrange(lDataplot, rDataplot,
                 labels = c("L", "R"), nCol = 2, common.legend = TRUE,
                 legend = "top")
taskFig


# summary stats for each task
library(Rmisc)
meanR_AE <- summarySE(data=rData, measurevar = "abs_err", groupvars = c("SUB", "TYPE"))

# just CL norm and CL fix
lData_NF <- lData[lData$TASK %in% c("cln","clf"), ]
rData_NF <- rData[rData$TASK %in% c("cln","clf"), ]

#plotting
ggplot(rData_NF, aes(x=SUB, y=abs_err, colour = TASK)) +
  geom_boxplot(outlier.alpha=0) +
  geom_jitter(position=position_dodge(width=.8), size=4, alpha=.5) +
  ylim(0,300) + geom_hline(yintercept=55, linetype="dotted") +
  theme_bw() + theme(legend.position = "bottom") -> rDataplot 
rDataplot
#same for left
ggplot(lData_NF, aes(x=SUB, y=abs_err, colour = TASK)) +
  geom_boxplot(outlier.alpha=0) +
  geom_jitter(position=position_dodge(width=.8), size=4, alpha=.5) +
  ylim(0,300) + geom_hline(yintercept=55, linetype="dotted") +
  theme_bw() + theme(legend.position = "bottom") -> lDataplot 
lDataplot

library(ggpubr)
taskFig_NF <- ggarrange(lDataplot, rDataplot,
                     labels = c("L", "R"), nCol = 2, common.legend = TRUE,
                     legend = "top")
taskFig_NF
