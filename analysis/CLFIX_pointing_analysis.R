library(readr)
library(ggplot2)
library(Rmisc)

#MANUALLY SET WORKING DIRECTORY TO DATA DIRECTORY AND LIST *.DAT CONTENTS
dataPath ="/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/DMT2019_analysis/pointing"
setwd(dataPath) #setting the path

#make complete list of eyelink RESULTS files (trial information)
resfiles=list.files(path=dataPath, full.names=TRUE, recursive = FALSE, include.dirs = TRUE, pattern = "*.csv")

#make empty dataframe to stack up trial information
CLF <- read.csv(text = "avg_rt,count_committed,count_valid_trial,datetime,eye_move,fix_x,land_x,land_y,reach_duration,response_time_target_mouse_response,sanity_check,subject_nr,subject_parity,targ_x,targ_y,target_onset,targets_file,time_touch_offset,touch_x")


#stack them all up
for(x in resfiles){
  tmp <- read.csv(x)
  #tmp$subject_nr <- as.numeric(substr(x))
  CLF <- rbind(CLF, tmp)
}

CLF <- CLF[CLF$targ_x > 0, c(5,7,8,12,14,15)]
CLF$GRP <- "HC"
CLF[CLF$subject_nr %in% c(801,804), "GRP"] <- "MCI"
CLF[CLF$subject_nr %in% c(802,803), "GRP"] <- "AD"
CLF$SUB <- factor(CLF$subject_nr)

CLF$x_Err <- CLF$land_x-CLF$targ_x #x and y error
CLF$y_Err <- CLF$land_y-CLF$targ_y
#converting to mm
pix_permm = 6.18249
CLF$x_Err = CLF$x_Err/pix_permm
CLF$y_Err = CLF$y_Err/pix_permm

CLF[is.na(CLF$eye_move), "eye_move"] <- 0 #any 'na' on eye-move = 0

CLF$ECC <- factor(cut(CLF$targ_x, 3, label=c("N", "M", "F")))
CLF <- CLF[CLF$eye_move == 0, ]

CLF$ABS_ERR <- sqrt(CLF$x_Err^2 + CLF$y_Err^2) #absolute error calc

#plotting individual participant data
ggplot(CLF, aes(x=SUB, y=ABS_ERR, colour = ECC)) +
  geom_boxplot(outlier.alpha=0) +
  geom_jitter(position=position_dodge(width=.8), size=4, alpha=.5) +
  ylim(0,50) + geom_hline(yintercept=10, linetype="dotted") +
  labs(x = 'Group', y = 'Absolute error (mm)') + 
  theme_bw() + theme(legend.position = "bottom") -> ppPlot
ppPlot

ggplot(CLF, aes(x=SUB, y=ABS_ERR, colour=ECC)) +
  geom_boxplot(outlier.alpha=0) +
  geom_jitter(position=position_dodge(width=.8), size=4, alpha=.5) +
  ylim(0,50) + geom_hline(yintercept=10, linetype="dotted")+
  theme_bw() + 


mean_AE <- summarySE(data=CLF, measurevar = "ABS_ERR", groupvars = c("SUB", "GRP"))

ggplot(mean_AE, aes(x = SUB, y=ABS_ERR, colour=GRP)) + geom_point(size=5, alpha=.5) +
  geom_errorbar(aes(ymin=ABS_ERR-ci, ymax=ABS_ERR+ci)) + ylim(c(50,175))-> ALL

grid.arrange(ALL, NOTNEAR, ncol=2)

mean_CE <- summarySE(data=CLF, measurevar = "x_Err", groupvars = c("SUB", "GRP"))

ggplot(mean_CE, aes(x = SUB, y=x_Err, colour=GRP)) + geom_point(size=5, alpha=.5) +
  geom_errorbar(aes(ymin=x_Err-ci, ymax=x_Err+ci))

lm_fits <- read.csv(text="SUB,rsq,a,b")

row=1

for(SUB in levels(CLF$SUB)){
  tmp <- CLF[CLF$SUB==SUB, ]
  model <- lm(ABS_ERR~ECC, data=tmp)
  lm_fits[row, 1] <- SUB
  lm_fits[row, 3:4] <- coefficients(model)[1:2]
  lm_fits[row, 2] <- summary(model)$r.squared
  row=row+1
}


#summary statistics
#all eccentricities
med <- aggregate(ABS_ERR~SUB*ECC, median, data=CLF) #median of eccentricity
mean_med <- aggregate(ABS_ERR~SUB, mean, data=med) #mean of all

ggplot(mean_med, aes(x=SUB, y=ABS_ERR))+ geom_point(size=5)

#far targets only
med_far <- aggregate(ABS_ERR~SUB*ECC, median, data=CLF[CLF$targ_x > 300, ]) #median of eccentricity
mean_med_far <- aggregate(ABS_ERR~SUB, mean, data=med) #mean of all

ggplot(mean_med_far, aes(x=SUB, y=ABS_ERR))+ geom_point(size=5)


