###########################
#ANALYSIS STAGE 1; reload AiA

#upload full dataframe
AiA <- read.csv("AiA_data.csv")

#define factors
AiA$ECC <- factor(AiA$ECC)
AiA$target <- factor(AiA$target)
AiA$TASK <- factor(AiA$TASK)
AiA$PPT <- factor(AiA$PPT)
AiA$LOOP <- factor(AiA$LOOP)

#DATA EXCLUSIONS
#exclude broken fixations
AiA <- AiA[AiA$RESP != -1, ]

#exclude movements ending in final frame (where this information is available)
AiA$FILT <- AiA$RT + AiA$MT > 2980
AiA[is.na(AiA$FILT), "FILT"] <- FALSE
AiA <- AiA[AiA$FILT == FALSE, ]

#exclude 'anticipations' (at RT < 1160, as 100 ms after target onset)
AiA$FILT <- AiA$RT < 1160
AiA[is.na(AiA$FILT), "FILT"] <- FALSE
AiA <- AiA[AiA$FILT == FALSE, ]

#remove filter column
AiA$FILT <- NULL

#visual screen for trial alignment
ggplot(AiA, aes(x=tX, y=rX, colour=TASK)) + geom_point(size=2, alpha =.3) + facet_wrap(~PPT)

#exclude participant 11 on basis of trial misalignment
AiA <- AiA[AiA$PPT != 11, ]

#check dual task performance
aggregate(DUAL~PPT, mean, data=AiA[AiA$TASK=="DUAL", ])

#EXCLUDING PARTICIPANTS 6 & 10 BASED ON DUAL TASK PERFORMANCE (accuracy < 0.5)
AiA<- AiA %>% filter(PPT != 10)
AiA<- AiA %>% filter(PPT != 6)

#relevel participant factor
AiA$PPT <- factor(AiA$PPT)


###########################
#ANALYSIS STAGE 2; plots and ANOVAs

#ANGULAR ERROR
#FULL PLOT
med_ANG <- aggregate(ANG_ERR~TASK*LOOP*ECC*PPT, mean, data=AiA)
plot_ANGERR <- summarySE(data=med_ANG, measurevar = "ANG_ERR", groupvars = c("TASK", "LOOP", "ECC"))

ggplot(plot_ANGERR, aes(x=ECC, y=ANG_ERR, colour=TASK, group=TASK)) +
  geom_point(size=5, alpha=.5, position=position_dodge(width=.3)) +
  geom_errorbar(aes(ymin=ANG_ERR-ci, ymax=ANG_ERR+ci), width=.4, position=position_dodge(width=.3)) +
  geom_line(position=position_dodge(width=.3)) +
  geom_hline(yintercept=0, linetype="dotted") +
  facet_wrap(~LOOP) +
  theme_bw()

#DIFFERENTIAL PLOT
wide_ANG <- dcast(PPT*LOOP*ECC~TASK, value.var = "ANG_ERR", data=med_ANG)
wide_ANG$DIFF <- wide_ANG$DUAL-wide_ANG$SINGLE
plot_ANG <- summarySE(data=wide_ANG, measurevar = "DIFF", groupvars = c("LOOP", "ECC"))

ggplot(plot_ANG, aes(x=ECC, y=DIFF, colour=LOOP, group=LOOP)) +
  geom_point(size=5, alpha=.5) +
  geom_errorbar(aes(ymin=DIFF-ci, ymax=DIFF+ci), width=.4) +
  geom_line() +
  geom_hline(yintercept=0, linetype="dotted") +
  ylim(c(-3,2)) +
  facet_wrap(~LOOP) +
  theme_bw()

#FULL ANOVA
ANG_ANOVA <- ezANOVA(
  data = med_ANG
  , dv = .(ANG_ERR)
  , wid = .(PPT)
  , within = .(TASK, LOOP, ECC)
  , type = 3
)

print(ANG_ANOVA)

#CLOSED-LOOP
ANG_CL <- ezANOVA(
  data = med_ANG[med_ANG$LOOP=="CL", ]
  , dv = .(ANG_ERR)
  , wid = .(PPT)
  , within = .(TASK, ECC)
  , type = 3
)

print(ANG_CL)


#OPEN-LOOP
rm_OL <- ezANOVA(
  data = med_ANG[med_ANG$LOOP=="OL", ]
  , dv = .(ANG_ERR)
  , wid = .(PPT)
  , within = .(TASK, ECC)
  , type = 3
)

print(rm_OL)



#AMPLITUDE ERROR
#FULL PLOT
med_AMP <- aggregate(AMP_ERR~TASK*LOOP*ECC*PPT, median, data=AiA)
plot_AMPERR <- summarySE(data=med_AMP, measurevar = "AMP_ERR", groupvars = c("TASK", "LOOP", "ECC"))

ggplot(plot_AMPERR, aes(x=ECC, y=AMP_ERR, colour=TASK, group=TASK)) +
  geom_point(size=5, alpha=.5, position=position_dodge(width=.3)) +
  geom_errorbar(aes(ymin=AMP_ERR-ci, ymax=AMP_ERR+ci), width=.4, position=position_dodge(width=.3)) +
  geom_line(position=position_dodge(width=.3)) +
  geom_hline(yintercept=0, linetype="dotted") +
  facet_wrap(~LOOP) +
  theme_bw()

#DIFFERENTIAL PLOT
wide_AMP <- dcast(PPT*LOOP*ECC~TASK, value.var = "AMP_ERR", data=med_AMP)
wide_AMP$DIFF <- wide_AMP$DUAL-wide_AMP$SINGLE
plot_AMP <- summarySE(data=wide_AMP, measurevar = "DIFF", groupvars = c("LOOP", "ECC"))

ggplot(plot_AMP, aes(x=ECC, y=DIFF, colour=LOOP, group=LOOP)) +
  geom_point(size=5, alpha=.5) +
  geom_errorbar(aes(ymin=DIFF-ci, ymax=DIFF+ci), width=.4) +
  geom_line(position=position_dodge(width=.2)) +
  geom_hline(yintercept=0, linetype="dotted") +
  facet_wrap(~LOOP) +
  theme_bw()

#FULL ANOVA
AMP_ANOVA <- ezANOVA(
  data = med_AMP
  , dv = .(AMP_ERR)
  , wid = .(PPT)
  , within = .(TASK, LOOP, ECC)
  , type = 3
)

print(AMP_ANOVA)

#CLOSED-LOOP
AMP_CL <- ezANOVA(
  data = med_AMP[med_AMP$LOOP=="CL", ]
  , dv = .(AMP_ERR)
  , wid = .(PPT)
  , within = .(TASK, ECC)
  , type = 3
)

print(AMP_CL)


#OPEN-LOOP
AMP_OL <- ezANOVA(
  data = med_AMP[med_AMP$LOOP=="OL", ]
  , dv = .(AMP_ERR)
  , wid = .(PPT)
  , within = .(TASK, ECC)
  , type = 3
)

print(AMP_OL)