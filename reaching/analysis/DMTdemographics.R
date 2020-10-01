## DMT 2019/20 REACHING IN AD ##
###### DEMOGRAPHICS DESCRIPTIVES ######
library(readr)
library(reshape2)
library(Rmisc)
library(ez)

##### PATHS + DATA #####
dataPath <- '/Users/Alex/Documents/DMT/data/'
setwd(dataPath)
patientDemo <- read.csv('patient_demographics.csv')
controlDemo <- read.csv('control_demographics.csv')

## adding site based on participant number
controlDemo$SITE <- factor(substr(controlDemo$subject_nr, 1, 1))
controlDemo$SITE <- factor(controlDemo$SITE, labels = c('UOE','UEA'))
patientDemo$SITE <- factor(substr(patientDemo$subject_nr, 1, 1))
patientDemo$SITE <- factor(patientDemo$SITE, labels = c('UOE','UEA'))

# remove 311
controlDemo <- controlDemo[controlDemo$subject_nr != 311 ,]

# get patient demos that match control and combine (for comparisons)
compDemo <- patientDemo[, c(1:6,18)]
compDemo <- rbind(controlDemo, compDemo)

##### SUMMARY STATS #####
## control + patient
# age
AGE <- summarySE(data = compDemo, measurevar = 'age', groupvars = c('diagnosis','SITE'))
# ANOVA
AGE_ANOVA <- ezANOVA(
  data = compDemo
  , dv = .(age)
  , wid = .(subject_nr)
  , between = .(SITE, diagnosis)
  , type = 3
)
print(AGE_ANOVA)

#pair-wise t-test
AGEttest <- pairwise.t.test(compDemo$age, compDemo$diagnosis, p.adj = 'bonf')
print(AGEttest)

# education
compDemo$education <- as.numeric(as.character(compDemo$education))
ED <- summarySE(data = compDemo, measurevar = 'education', groupvars = c('diagnosis','SITE'),
                na.rm = TRUE)
# ANOVA
# remove NAs
edAnova <- na.omit(compDemo)
ED_ANOVA <- ezANOVA(
  data = edAnova
  , dv = .(education)
  , wid = .(subject_nr)
  , between = .(SITE, diagnosis)
  , type = 3
)
print(ED_ANOVA)

## patient only 
# ACE
ACE <- summarySE(data = patientDemo, measurevar = 'ACEall', groupvars = c('diagnosis','SITE'))
# time since diagnosis
DAYS <- summarySE(data = patientDemo, measurevar = 'Days.since.diagnosis', 
                  groupvars = c('diagnosis','SITE'), na.rm = TRUE)
WEEKS <- summarySE(data = patientDemo, measurevar = 'Weeks.since.diagnosis', 
                   groupvars = c('diagnosis','SITE'), na.rm = TRUE)

