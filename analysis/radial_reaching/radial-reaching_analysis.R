library(readr)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(Rmisc)
library(ez)
library(psychReport)
library(singcar)

#on mac
#anaPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/analysis/radial_reaching'
#UEAPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/norwich_movement_data'
#dataPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/data'

# on desktop mac
anaPath <- '/Users/Alex/Documents/DMT/analysis/radial_reaching'
dataPath <- '/Users/Alex/Documents/DMT/data'
UEAPath <- '/Users/Alex/Documents/DMT/norwich_movement_data/'

#on pc
#dataPath <- 'S:/groups/DMT/data'
#anaPath <- 'S:/groups/DMT/analysis/radial_reaching'
setwd(anaPath) #for Edinburgh data 
resUOE <- read.csv('radial-reaching_compiled.csv')
setwd(UEAPath) #for Norwich data
resUEA <- read.csv('radial-reaching_compiled.csv') 

##### DATA ORGANISE #####
# UOE data, getting rid of 'deg' values - just use mm
resUOE <- resUOE[, c(1:7,14,15,8:13,16:19,22,24,25)]
# UEA data getting rid of not needed values
resUEA <- resUEA[, c(1:9,15:18,21,22,24:30)]
resUEA <- resUEA[resUEA$PPT != '311' ,]
# change group labelling in UEA (3 = control, 4 = patient) to match UOE (1 = control, 2 = patient)
resUEA$GRP <- factor(resUEA$GRP, labels = c('1','2'))
# change position labelling here too - needs to match UOE (100,200,300,400mm)
# split into two data-frames (left/right sides) and label, the bind again
UEAright <- resUEA[resUEA$SIDE == 'Right' ,]
UEAleft <- resUEA[resUEA$SIDE == 'Left' ,]
# for right 1= closest, 4= furthest away
UEAright$POSITION <- factor(UEAright$POSITION, labels = c('100','200','300','400'))
# for left 1=furthest, 4= nearest
UEAleft$POSITION <- factor(UEAleft$POSITION, labels = c('-400','-300','-200','-100'))
# bind back together!
resUEA <- rbind(UEAright, UEAleft)

# now let's bind! & order by participant
res <- rbind(resUOE, resUEA)
res <- res[order(res$PPT),]

# changing levelsfor plotting
res$VIEW <- factor(res$VIEW) #changing so only 2 levels recorded
res$GRP <- factor(res$GRP)
levels(res$GRP) <- c('Control', 'Patient') #changing group name from 1 = control, 2 = AD
res$diagnosis <- factor(res$diagnosis)
res$POSITION <- factor(res$POSITION)
# capitalising names of demographics to match rest
colnames(res)[colnames(res) == 'gender'] <- 'GENDER'
colnames(res)[colnames(res) == 'age'] <- 'AGE'
colnames(res)[colnames(res) == 'hand'] <- 'HAND'
colnames(res)[colnames(res) == 'education'] <- 'ED'
colnames(res)[colnames(res) == 'diagnosis'] <- 'DIAGNOSIS'

## getting dominant + non-dominant sides, for analysis
res$HAND <- factor(res$HAND, labels = c('Left','Right'))
# if hand = side, dominant; else non dominant
res$DOM <- as.numeric(res$SIDE == res$HAND) #1 = dominant, 0 = non-dominant
res$DOM <- factor(res$DOM, labels= c('ND','D'))
# change order so important var up-front
res <- res[, c(1,6,21:23,2:5,7:20)]
# save Edinburgh & UEA compiled data
setwd(anaPath)
write.csv(res, 'all_radial-reaching_compiled.csv', row.names = FALSE)

# summary data
# extracting data from furthest two target locs
res_periph <- subset(res, res$POSITION == -200 | res$POSITION == -300 | 
                         res$POSITION == 300 | res$POSITION == 200)
res_medians <- aggregate(AE ~ PPT*POSITION*VIEW*SIDE*DOM*DIAGNOSIS*GRP*SITE*AGE*ED, 
                         mean, data = res_periph)
colnames(res_medians)[colnames(res_medians)=='AE'] <- 'AEmed'

# removing free + peripheral trails at 100mm left for 101, error
res_medians <- res_medians[order(res_medians$PPT) ,]
res_medians <- res_medians[!(res_medians$PPT == '101' & res_medians$AEmed > 50) ,]

res_means <- aggregate(AEmed ~ PPT*VIEW*SIDE*DOM*GRP*DIAGNOSIS*SITE*AGE*ED, 
                       mean, data = res_medians)
colnames(res_means)[colnames(res_means)=='AEmed'] <- 'AEmean'

# casting by task
PMIdata <- dcast(res_means, PPT+GRP+DIAGNOSIS+DOM+SIDE+SITE+AGE+ED ~ VIEW)
PMIdata$PMI <- PMIdata$Peripheral - PMIdata$Free

## summary data
meanPMI_side <- summarySE(PMIdata, measurevar = 'PMI', groupvar = c('DIAGNOSIS', 'DOM'),
                          na.rm = TRUE)
meanPMI_all <- summarySE(PMIdata, measurevar = 'PMI', groupvar = c('DIAGNOSIS'),
                         na.rm = TRUE)

##### ANOVA - differences between sites #####
siteANOVA <- na.omit(res_means)
siteANOVA <- siteANOVA[siteANOVA$PPT != 212 ,]
siteANOVA <- siteANOVA[siteANOVA$PPT != 310 , c(1:7,10)]
siteANOVA <- aggregate(AEmean ~ PPT+VIEW+DOM+DIAGNOSIS+SITE, mean, data = siteANOVA)

SITE_ANOVA <- ezANOVA(
  data = siteANOVA
  , dv = .(AEmean)
  , wid = .(PPT)
  , within = .(VIEW)
  , between =. (SITE)
  , type = 3
)

print(SITE_ANOVA)
# no difference across sites! Hooray!

ggplot(siteANOVA, aes(x = VIEW, y = AEmean, colour = DIAGNOSIS)) +
  geom_point(size =2, position = position_dodge(.3)) +
  facet_wrap(~SITE)

ggsave('site_comparison.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

# Eccentricity plots
# creating another data-frame with all position data, for plotting eccentricity info
res_medians_all <- aggregate(AE ~ PPT*POSITION*VIEW*SIDE*DOM*DIAGNOSIS*GRP*SITE*AGE*ED, 
                             mean, data = res)
# removing free + peripheral trails at 100mm left for 101, error
res_medians_all <- res_medians_all[order(res_medians_all$PPT) ,]
res_medians_all <- res_medians_all[!(res_medians_all$PPT == '101' & res_medians_all$AE > 50) ,]

# controls
#re-ordering position data to fit correct order
res_medians_all$POSITION <- as.numeric(as.character(res_medians_all$POSITION))
res_medians_all$POSITION <- factor(res_medians_all$POSITION)
meds_control <- res_medians_all[res_medians_all$GRP == 'Control' ,]

ggplot(meds_control, aes(x = POSITION, y = AE, colour = DOM)) +
  geom_point(shape = 1, size = 2) +
  geom_line(aes(group = PPT), size = 0.5, alpha = .5) +
  facet_wrap(~VIEW) + 
  labs(title = 'Control', x = 'Eccentricity (mm)', 
       y = 'Mean AE (mm)', element_text(size = 12)) +
  scale_colour_manual(values = c('black', 'grey50')) +
  theme_bw() + theme(legend.position = 'none', text = element_text(size = 10),
                     strip.text.x = element_text(size = 10)) 


ggsave('control_ecc.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

# patients - MCI
meds_MCI <- res_medians_all[res_medians_all$DIAGNOSIS == 'MCI' ,]

ggplot(meds_MCI, aes(x = POSITION, y = AE, colour = SIDE)) +
  geom_point(shape = 1, size = 2) +
  geom_line(aes(group = PPT), size = 0.5, alpha = .5) +
  facet_wrap(~VIEW) +
  labs(title = 'MCI', x = 'Eccentricity (mm)', 
       y = 'Mean AE (mm)', element_text(size = 12)) +
  scale_colour_manual(values = c('black', 'grey50')) +
  theme_bw() + theme(legend.position = 'none', text = element_text(size = 10),
                     strip.text.x = element_text(size = 10)) 

ggsave('MCI_ecc.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

# patients - AD
meds_AD <- res_medians_all[res_medians_all$DIAGNOSIS == 'AD' ,]

ggplot(meds_AD, aes(x = POSITION, y = AE, colour = SIDE)) +
  geom_point(shape = 1, size = 2) +
  geom_line(aes(group = PPT), size = 0.5, alpha = .5) +
  facet_wrap(~VIEW) + 
  labs(title = 'Alzheimers', x = 'Eccentricity (deg)', 
       y = 'Mean AE (mm)', element_text(size = 12)) +
  scale_colour_manual(values = c('black', 'grey50')) +
  theme_bw() + theme(legend.position = 'none', text = element_text(size = 10),
                     strip.text.x = element_text(size = 10)) 

ggsave('AD_ecc.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, path = anaPath)

# correlate PMI with possible co-variates, age + education
agecov <- cor.test(PMIdata$AGE, PMIdata$PMI, method = 'pearson')
ggscatter(PMIdata, x = "AGE", y = "PMI", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Age", ylab = "PMI (mm)")

# correlate PMI with years of education
PMIdata$ED <- as.numeric(PMIdata$ED)
edcov <- cor.test(PMIdata$ED, PMIdata$PMI, method = 'pearson')
ggscatter(PMIdata, x = "ED", y = "PMI", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Years of Education", ylab = "PMI (mm)")

######## OUTLIER REMOVAL, filtered PMI #########
# two different sites and set-ups, use data from control group in each site to remove outliers
# need to run same removal for outliers in each site
## UOE FIRST ##
controlUOE <- PMIdata[PMIdata$PPT < 200, ]

# median values for each side
tmp <- aggregate(PMI ~ GRP, median, data = controlUOE)
names(tmp)[2] <- 'med'
controlUOE <- merge(tmp, controlUOE)

# calculating MAD for each (absolute value)
controlUOE$AD <- abs(controlUOE$PMI - controlUOE$med)
tmp <- aggregate(AD ~ GRP, median, data=controlUOE)
names(tmp)[2] <- 'MAD'
controlUOE <- merge(controlUOE, tmp)

# adjusted z-score from these values
controlUOE$az <- (controlUOE$PMI - controlUOE$med)/(controlUOE$MAD * 1.4826)
controlUOE$z <- scale(controlUOE$PMI)
controlUOE$PPT <- factor(controlUOE$PPT)

plot_name = 'adjustedZ_UOE.png'
ggplot(controlUOE, aes(x = GRP, y = az, colour = PPT)) +
  geom_point(size = 3, position = position_dodge(.1)) +  
  stat_summary(aes(y = az, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', group = 1) + 
  labs(x = '', y = 'Adjusted z-score', element_text(size = 13)) +
  theme_bw() + theme(legend.position = "right", legend.title = element_blank())

ggsave(plot_name, plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, width = 5, height = 6.5, path = anaPath)

# find controls with az > 2.5 - need to remove entire control, not just side
XCLUDE_UOE <- controlUOE[controlUOE$az > 2.5, ]

## UEA SECOND ##
controlUEA <- PMIdata[PMIdata$SITE == 'UEA' ,]
controlUEA <- controlUEA[controlUEA$PPT < 400 ,]

# median values for each side
tmp <- aggregate(PMI ~ GRP, median, data = controlUEA)
names(tmp)[2] <- 'med'
controlUEA <- merge(tmp, controlUEA)

# calculating MAD for each (absolute value)
controlUEA$AD <- abs(controlUEA$PMI - controlUEA$med)
tmp <- aggregate(AD ~ GRP, median, data=controlUEA)
names(tmp)[2] <- 'MAD'
controlUEA <- merge(controlUEA, tmp)

# adjusted z-score from these values
controlUEA$az <- (controlUEA$PMI - controlUEA$med)/(controlUEA$MAD * 1.4826)
controlUEA$z <- scale(controlUEA$PMI)
controlUEA$PPT <- factor(controlUEA$PPT)

plot_name = 'adjustedZ_UEA.png'
ggplot(controlUEA, aes(x = GRP, y = az, colour = PPT)) +
  geom_point(size = 3, position = position_dodge(.1)) +  
  stat_summary(aes(y = az, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', group = 1) + 
  labs(x = '', y = 'Adjusted z-score', element_text(size = 13)) +
  theme_bw() + theme(legend.position = "right", legend.title = element_blank())

ggsave(plot_name, plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, width = 5, height = 6.5, path = anaPath)

# find controls with az > 2.5 - need to remove entire control, not just side
XCLUDE_UEA <- controlUEA[controlUEA$az > 2.5, ]

### combinging outliers and removing
setwd(anaPath)
XCLUDE <- rbind(XCLUDE_UOE,XCLUDE_UEA)
write.csv(XCLUDE, 'radial-reaching_outliers.csv', row.names = FALSE)
# creating data-frame with control data removed
PMIfilt <- PMIdata[!(PMIdata$PPT %in% XCLUDE$PPT), ]
# get rid of NA values
PMIfilt <- na.omit(PMIfilt)

res_mediansF <- res_medians[!(res_medians$PPT %in% XCLUDE$PPT), ]
res_medians_allF <- res_medians_all[!(res_medians_all$PPT %in% XCLUDE$PPT), ]
res_meansF <- res_means[!(res_means$PPT %in% XCLUDE$PPT), ]
# saving filered data
setwd(anaPath)
write.csv(PMIfilt, 'radialPMI-filtered.csv', row.names = FALSE)
write.csv(res_mediansF, 'radial-medians_filtered.csv', row.names = FALSE)
write.csv(res_meansF, 'radial-means_filtered.csv', row.names = FALSE)
write.csv(res_medians_allF, 'radial-medians-all_filtered.csv', row.names = FALSE)

meanPMI_side <- summarySE(PMIfilt, measurevar = 'PMI', groupvar = c('DIAGNOSIS', 'DOM'),
                          na.rm = TRUE)
meanPMI_all <- summarySE(PMIfilt, measurevar = 'PMI', groupvar = c('DIAGNOSIS'),
                         na.rm = TRUE)

######## CASE CONTROL ANALYSIS ########
# create data-frames (controls and patients) with key information
td_summary <- summarySE(data = PMIfilt, measurevar = 'PMI', 
                        groupvars = c('DIAGNOSIS','DOM','SITE'), na.rm = TRUE)
#patient data for test of deficit
td_patient <- PMIfilt[, c(1,3,4,6,11)]
td_patient <- td_patient[td_patient$DIAGNOSIS != 'HC' ,] #removing controls
# next, run case-control independently for each site - different set ups

## UOE CASE CONTROL ##
# isolating control data for analysis
tdUOE_control <- td_summary[td_summary$DIAGNOSIS == 'HC' 
                            & td_summary$SITE == 'UOE' ,]
tdUOE_patient <- td_patient[td_patient$SITE == 'UOE' ,]
tdUOE_patient <- dcast(tdUOE_patient, PPT+DIAGNOSIS ~ DOM)
# NA values  = -1, so t-value is negative and can remove later
tdUOE_patient[is.na(tdUOE_patient$D), "D"] <- -1 

# time for test of deficit! Calling on 'singcar' package developed by Jonathan Rittmo
# using Crawford's (1998) test of deficit
td_dom <- read.csv(text = 'PMI,TSTAT,PVALUE,PROP,ZCC,CI,LCI-T,HCI-T,LCI-P,HCI-P,DOM,PPT,DIAGNOSIS')
td_ndom <- read.csv(text = 'PMI,TSTAT,PVALUE,PROP,ZCC,CI,LCI-T,HCI-T,LCI-P,HCI-P,DOM,PPT,DIAGNOSIS')
for (l in 1:length(tdUOE_patient$PPT)){
  #left data first
  leftres <- TD(tdUOE_patient$ND[l], tdUOE_control$PMI[1], tdUOE_control$sd[1], 24, 
                alternative = 'greater', na.rm = FALSE)
  diff <- t(leftres$estimate)
  ltmp <- data.frame(tdUOE_patient$ND[l], leftres$statistic, leftres$p.value, 
                     diff[1], diff[2], t(leftres$interval), 'ND', check.names = FALSE) 
  ltmp$PPT <- tdUOE_patient$PPT[l]
  ltmp$DIAGNOSIS <- tdUOE_patient$DIAGNOSIS[l]
  td_ndom <- rbind(td_ndom, ltmp)
  #then right data
  rightres <- TD(tdUOE_patient$D[l], tdUOE_control$PMI[2], tdUOE_control$sd[2], 24, 
                 alternative = 'greater', na.rm = FALSE)
  diff <- t(rightres$estimate)
  rtmp <- data.frame(tdUOE_patient$D[l], rightres$statistic, rightres$p.value, 
                     diff[1], diff[2], t(rightres$interval), 'D', check.names = FALSE) 
  rtmp$PPT <- tdUOE_patient$PPT[l]
  rtmp$DIAGNOSIS <- tdUOE_patient$DIAGNOSIS[l]
  td_dom <- rbind(td_dom, rtmp)
}

# merging and renaming data-frames
tdUOE_results <- read.csv(text = 'PMI,TSTAT,PVALUE,ZCC,PROP,CI,LCI-T,HCI-T,LCI-P,HCI-P,DOM,PPT,DIAGNOSIS')
# changing names of td data-frames to match td-res
names(td_ndom) <- names(tdUOE_results)
names(td_dom) <- names(tdUOE_results)
tdUOE_results <- rbind(tdUOE_results, td_ndom, td_dom)
# remove original NA values (PMI = -1), where patients had no data
tdUOE_results <- tdUOE_results[tdUOE_results$PMI > -0.9 ,]
tdUOE_results$PPT <- factor(tdUOE_results$PPT)
# add site back to results
tdUOE_results$SITE <- 'UOE'

# identifying patients with significant deficit
tdUOE_results$DEFICIT <- tdUOE_results$PVALUE < 0.025
# identifying patients with possible deficit
tdUOE_results$BL <- tdUOE_results$PVALUE < 0.05
# convert logicals into numerics, useful for binomial later
tdUOE_results$DEFICIT <- as.numeric(tdUOE_results$DEFICIT)
tdUOE_results$BL <- as.numeric(tdUOE_results$BL)

## UEA CASE CONTROL ##
# isolating control data for analysis
tdUEA_control <- td_summary[td_summary$DIAGNOSIS == 'HC' 
                            & td_summary$SITE == 'UEA' ,]
tdUEA_patient <- td_patient[td_patient$SITE == 'UEA' ,]
tdUEA_patient <- dcast(tdUEA_patient, PPT+DIAGNOSIS ~ DOM)

# time for test of deficit! Calling on 'singcar' package developed by Jonathan Rittmo
# using Crawford's (1998) test of deficit
td_dom <- read.csv(text = 'PMI,TSTAT,PVALUE,PROP,ZCC,CI,LCI-T,HCI-T,LCI-P,HCI-P,DOM,PPT,DIAGNOSIS')
td_ndom <- read.csv(text = 'PMI,TSTAT,PVALUE,PROP,ZCC,CI,LCI-T,HCI-T,LCI-P,HCI-P,DOM,PPT,DIAGNOSIS')
for (l in 1:length(tdUEA_patient$PPT)){
  #left data first
  leftres <- TD(tdUEA_patient$ND[l], tdUEA_control$PMI[1], tdUEA_control$sd[1], 24, 
                alternative = 'greater', na.rm = FALSE)
  diff <- t(leftres$estimate)
  ltmp <- data.frame(tdUEA_patient$ND[l], leftres$statistic, leftres$p.value, 
                     diff[1], diff[2], t(leftres$interval), 'ND', check.names = FALSE) 
  ltmp$PPT <- tdUEA_patient$PPT[l]
  ltmp$DIAGNOSIS <- tdUEA_patient$DIAGNOSIS[l]
  td_ndom <- rbind(td_ndom, ltmp)
  #then right data
  rightres <- TD(tdUEA_patient$D[l], tdUEA_control$PMI[2], tdUEA_control$sd[2], 24, 
                 alternative = 'greater', na.rm = FALSE)
  diff <- t(rightres$estimate)
  rtmp <- data.frame(tdUEA_patient$D[l], rightres$statistic, rightres$p.value, 
                     diff[1], diff[2], t(rightres$interval), 'D', check.names = FALSE) 
  rtmp$PPT <- tdUEA_patient$PPT[l]
  rtmp$DIAGNOSIS <- tdUEA_patient$DIAGNOSIS[l]
  td_dom <- rbind(td_dom, rtmp)
}

# merging and renaming data-frames
tdUEA_results <- read.csv(text = 'PMI,TSTAT,PVALUE,ZCC,PROP,CI,LCI-T,HCI-T,LCI-P,HCI-P,DOM,PPT,DIAGNOSIS')
# changing names of td data-frames to match td-res
names(td_ndom) <- names(tdUEA_results)
names(td_dom) <- names(tdUEA_results)
tdUEA_results <- rbind(tdUEA_results, td_ndom, td_dom)
# remove original NA values (PMI = -1), where patients had no data
tdUEA_results <- tdUEA_results[tdUEA_results$PMI > -0.9 ,]
tdUEA_results$PPT <- factor(tdUEA_results$PPT)
# add site back to results
tdUEA_results$SITE <- 'UEA'

# identifying patients with significant deficit
tdUEA_results$DEFICIT <- tdUEA_results$PVALUE < 0.025
# identifying patients with possible deficit
tdUEA_results$BL <- tdUEA_results$PVALUE < 0.05
# convert logicals into numerics, useful for binomial later
tdUEA_results$DEFICIT <- as.numeric(tdUEA_results$DEFICIT)
tdUEA_results$BL <- as.numeric(tdUEA_results$BL)

## merge UOE and UEA data-frames
td_results <- rbind(tdUOE_results, tdUEA_results)
#save
write.csv(td_results, 'radial-reaching_case-control.csv', row.names = FALSE)

##### BINOMIAL STATS #####
# calculating the likelihood that inflated PMI occurs at above chance in each patient group
# with filtered data. Deficit and borderline deficit
pval <- .05
pval_Bl <- .1

## take any p < .025 for either side, re-run binomial on this data
# extracting p-value for each side (left and right) using td_results
td_side <- dcast(td_results, PPT+DIAGNOSIS+SITE ~ DOM, value.var = 'PVALUE')
td_side[is.na(td_side$D), "D"] <- 10
# finding defcits and borderline cases
td_side$DEFICIT <- (td_side$ND < 0.025 | td_side$D < 0.025)
td_side$BL <- (td_side$ND < 0.05 | td_side$D < 0.05)
# changing to numerical value
td_side$DEFICIT <- as.numeric(td_side$DEFICIT)
td_side$BL <- as.numeric(td_side$BL)

# splitting into MCI and AD
MCIP <- td_side[td_side$DIAGNOSIS == 'MCI', ]
ADP <- td_side[td_side$DIAGNOSIS == 'AD', ]

# binomial test of deficit
binMCI <- binom.test(sum(MCIP$DEFICIT), length(MCIP$PPT), pval, alternative = 'greater')
print(binMCI)
binAD <- binom.test(sum(ADP$DEFICIT), length(ADP$PPT), pval, alternative = 'greater')
print(binAD)
binALL <- binom.test(sum(td_side$DEFICIT), length(td_side$PPT), pval, alternative = 'greater')
print(binALL)

# binomial test of borderline deficit
binMCI <- binom.test(sum(MCIP$BL), length(MCIP$PPT), pval_Bl, alternative = 'greater')
print(binMCI)
binAD <- binom.test(sum(ADP$BL), length(ADP$PPT), pval_Bl, alternative = 'greater')
print(binAD)
binALL <- binom.test(sum(td_side$BL), length(td_side$PPT), pval_Bl, alternative = 'greater')
print(binALL)

##### ANOVAS, group comparisons #####
# use PMIfilt data-frame to run between-subject ANOVA
# removing NA values for ANOVA - entire participant (not just side)
PMIanova <- PMIfilt[PMIfilt$PPT != 212, ]

# FULL ANOVA ON FILTERED DATA
FILT_ANOVA <- ezANOVA(
  data = PMIanova
  , dv = .(PMI)
  , wid = .(PPT)
  , within = .(DOM)
  , between = .(DIAGNOSIS)
  , between_covariates = .(SITE)
  , type = 3,
  return_aov = TRUE,
  detailed = TRUE
)

FILT_ANOVA$ANOVA
FILT_ANOVA$`Mauchly's Test for Sphericity`
FILT_ANOVA$`Sphericity Corrections`
aovPMI <- aovEffectSize(ezObj = FILT_ANOVA, effectSize = "pes")
aovDispTable(aovPMI)

#pair-wise t-test
PMIttest <- pairwise.t.test(PMIanova$PMI, PMIanova$DIAGNOSIS, p.adj = 'bonf')
print(PMIttest)

### ANOVA BY ECCENTRICITY ###
# converting factor back to numeric keeping values
tmp <- as.numeric(levels(res_medians_allF$POSITION))[res_medians_allF$POSITION]

res_medians_allF$ECC <- abs(tmp)
ECCanova <- res_medians_allF[res_medians_allF$ECC != 400 ,]

# removing participants with incomplete data-sets
ECCanova <- ECCanova[ECCanova$PPT != 212 & ECCanova$PPT != 101
                     & ECCanova$PPT != 310,]
ECCanova$ECC <- factor(ECCanova$ECC)

ECC_ANOVA <- ezANOVA(
  data = ECCanova
  , dv = .(AE)
  , wid = .(PPT)
  , within = .(DOM, VIEW, ECC)
  , between = .(DIAGNOSIS)
  , between_covariates = .(SITE)
  , type = 3,
  return_aov = TRUE,
  detailed = TRUE
)

ECC_ANOVA$ANOVA
ECC_ANOVA$`Mauchly's Test for Sphericity`
ECC_ANOVA$`Sphericity Corrections`
aovECC <- aovEffectSize(ezObj = ECC_ANOVA, effectSize = "pes")
aovDispTable(aovECC)

#pair-wise t-test
ECCttest <- pairwise.t.test(ECCanova$AE, ECCanova$DIAGNOSIS, p.adj = 'bonf')
print(ECCttest)

##### PLOTTING #####
## PLOT 1: eccentricity ##
ECCsummary <- summarySE(res_medians_allF, measurevar = 'AE', 
                        groupvar = c('DIAGNOSIS', 'SIDE', 'VIEW', 'POSITION'), na.rm = TRUE)
ECCsummary$DIAGNOSIS <- factor(ECCsummary$DIAGNOSIS, levels = c('HC','MCI','AD'))

ggplot(ECCsummary, aes(x = POSITION, y = AE, group = DIAGNOSIS, colour = DIAGNOSIS, 
                   shape = DIAGNOSIS)) +
  geom_point(size = 3, position = position_dodge(width = .4)) +
  geom_errorbar(aes(ymin=AE-ci, ymax=AE+ci), 
                width=.4, position = position_dodge(width = .4)) + 
  geom_line(aes(group = DIAGNOSIS), size = 0.7, position = position_dodge(width = .4)) +
  scale_color_manual(values = c('black','grey30','grey60')) +
  labs(x = 'Eccentricity (mm)', y = 'Radial reaching error (mm)') +
  facet_grid(~VIEW) + theme_classic() +
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12)
  )

ggsave('RADeccentricity-fig.png', plot = last_plot(), device = NULL, dpi = 300, 
       width = 7, height = 5, path = anaPath)

## PLOT 2: mean AE
AEsummary <- summarySE(res_meansF, measurevar = 'AEmean', 
                       groupvar = c('DIAGNOSIS', 'DOM', 'VIEW'), na.rm = TRUE)

AEsummary$DIAGNOSIS <- factor(AEsummary$DIAGNOSIS, levels = c('HC','MCI','AD'))

ggplot(AEsummary, aes(x = DOM, y = AEmean, colour = VIEW, group = DIAGNOSIS)) +
  geom_point(size = 3) +
  geom_line(aes(group = VIEW), size = 0.7) +
  geom_errorbar(aes(ymin=AEmean-ci, ymax=AEmean+ci), 
                width=.3) +
  scale_color_manual(values = c('black','grey60')) +
  labs(x = 'Side', y = 'Radial reaching error (mm)') + 
  facet_wrap(~DIAGNOSIS) +
  theme_classic() + theme(legend.position = 'bottom', 
                          legend.title = element_blank(),
                          axis.text = element_text(size = 10),
                          axis.title = element_text(size = 12),
                          strip.text = element_text(size = 12))

ggsave('RADmean-fig.png', plot = last_plot(), device = NULL, dpi = 300, 
       width = 4.5, height = 6, path = anaPath)

## PLOT 3: PMI ##
## PLOT 3: PMI ##
# make control data-frame
control_PMI <- subset(PMIfilt, PMIfilt$DIAGNOSIS == 'HC')
control_PMI$TSTAT <- 0
control_PMI$PVALUE <- 1
control_PMI$DEFICIT <- 0
control_PMI$BL <- 0
# get deficit data for patients
plot_PMI <- merge(PMIfilt, td_results, by = c('PPT','DOM','DIAGNOSIS','SITE'))
# include only relevant info
plot_PMI$PMI <- plot_PMI$PMI.x
plot_PMI <- plot_PMI[, c(1:10,24,13,14,22,23)]

# make plot data frame
plot_PMI <- rbind(control_PMI, plot_PMI)

plot_PMI$DOM <- factor(plot_PMI$DOM, levels = c('ND', 'D'))
plot_PMI$DIAGNOSIS <- factor(plot_PMI$DIAGNOSIS, levels = c('HC','MCI','AD'))
plot_PMI$PPT <- factor(plot_PMI$PPT)
# make column where deficit and BL cases are combined
# 1 = no deficit, 2 = borderline deficit, 3 = deficit
plot_PMI$DEFICITS <- as.numeric(plot_PMI$DEFICIT) + as.numeric(plot_PMI$BL)
plot_PMI$DEFICITS <- factor(plot_PMI$DEFICITS)

ggplot(plot_PMI, aes(x = DOM, y = PMI, colour = DIAGNOSIS, group = PPT, shape = DEFICITS)) + 
  geom_line(aes(group = PPT), alpha = .7, size = 0.7, position = position_dodge(.2)) +
  geom_point(size = 2.5, position = position_dodge(.2)) +
  stat_summary(aes(y = PMI, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1, size = 4, group = 1) +
  scale_color_manual(values = c('grey45','grey45','grey45')) +
  scale_shape_manual(values = c(1,18,16)) +
  facet_wrap(~DIAGNOSIS) +
  labs(x = 'Side', y = 'Lateral PMI (mm)') +
  theme_classic() + theme(legend.position = 'none', 
                          axis.title = element_text(size = 12),
                          axis.text = element_text(size = 10),
                          strip.text = element_text(size = 10) 
  ) -> pPMI
pPMI

## PLOT 4: average PMI across sides - combine with PLOT 2
PMIav_plot <- aggregate(PMI ~ PPT*DIAGNOSIS*AGE*ED*SITE, mean, data = plot_PMI)
PMIav_plot <- PMIav_plot[order(PMIav_plot$PPT),]
plot_PMI$DEFICITS <- as.numeric(plot_PMI$DEFICITS)
deficit <- aggregate(DEFICITS ~ PPT*DIAGNOSIS, max, data = plot_PMI)
deficit <- deficit[order(deficit$PPT),]
PMIav_plot$DEFICITS <- factor(deficit$DEFICITS)

ggplot(PMIav_plot, aes(x = DIAGNOSIS, y = PMI, colour = DIAGNOSIS, group = PPT, shape = DEFICITS)) + 
  geom_point(size = 2.5, position = position_dodge(.2)) +
  stat_summary(aes(y = PMI, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1, size = 4, group = 1) +
  scale_color_manual(values = c('grey45','grey45','grey45')) +
  scale_shape_manual(values = c(1,18,16)) +
  labs(x = 'Diagnosis', y = '') +
  theme_classic() + theme(legend.position = 'none', 
                          axis.title = element_text(size = 12),
                          axis.text = element_text(size = 10),
                          strip.text = element_text(size = 10) 
  ) -> avPMI
avPMI

PMIfig <- ggarrange(pPMI, avPMI,
                    ncol=2, nrow=1,
                    widths = c(1.5,1),
                    labels = c('a','b'),
                    hjust = -1)
PMIfig
ggsave('RADPMI-fig.png', plot = last_plot(), device = NULL, dpi = 300, 
       width = 8, height = 4, path = anaPath)

