library(readr)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(Rmisc)

#on mac
#anaPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/analysis/radial_reaching'
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
res_medians_all <- aggregate(AE ~ PPT*POSITION*VIEW*SIDE*DOM*DIAGNOSIS*GRP*SITE*AGE*ED, 
                         mean, data = res)
colnames(res_medians_all)[colnames(res_medians_all)=='AE'] <- 'AEmed'
# removing free + peripheral trails at 100mm left for 101, error
res_medians <- res_medians[order(res_medians$PPT) ,]
res_medians <- res_medians[!(res_medians$PPT == '101' & res_medians$AEmed > 50) ,]

res_means <- aggregate(AEmed ~ PPT*VIEW*SIDE*DOM*GRP*DIAGNOSIS*SITE*AGE*ED, 
                       mean, data = res_medians)
colnames(res_means)[colnames(res_means)=='AEmed'] <- 'AEmean'
write.csv(res_medians, 'radial-reaching_medians.csv', row.names = FALSE)
write.csv(res_means, 'radial-reaching_means.csv', row.names = FALSE)

# casting by task
PMIdata <- dcast(res_means, PPT+GRP+DIAGNOSIS+DOM+SIDE+SITE+AGE+ED ~ VIEW)
PMIdata$PMI <- PMIdata$Peripheral - PMIdata$Free
write.csv(PMIdata, 'radial-reaching_PMI.csv', row.names = FALSE)

## summary data
meanPMI_side <- summarySE(PMIdata, measurevar = 'PMI', groupvar = c('DIAGNOSIS', 'DOM'),
                          na.rm = TRUE)
meanPMI_all <- summarySE(PMIdata, measurevar = 'PMI', groupvar = c('DIAGNOSIS'),
                         na.rm = TRUE)

#### plot by eccentricity
# creating another data-frame with all position data, for plotting eccentricity info
res_medians_all <- aggregate(AE ~ PPT*POSITION*VIEW*SIDE*DOM*DIAGNOSIS*GRP*SITE*AGE*ED, 
                             mean, data = res)
# removing free + peripheral trails at 100mm left for 101, error
res_medians_all <- res_medians_all[order(res_medians_all$PPT) ,]
res_medians_all <- res_medians_all[!(res_medians_all$PPT == '101' & res_medians_all$AE > 50) ,]

# controls
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
#### UOE FIRST #####
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

### outlier removal ###
# find controls with az > 2.5 - need to remove entire control, not just side
XCLUDE_UOE <- controlUOE[controlUOE$az > 2.5, ]
write.csv(XCLUDE_UOE, 'radialoutliers_UOE.csv', row.names = FALSE)

#### UEA SECOND #####
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

### outlier removal ###
# find controls with az > 2.5 - need to remove entire control, not just side
XCLUDE_UEA <- controlUEA[controlUEA$az > 2.5, ]
write.csv(XCLUDE_UEA, 'radialoutliers_UEA.csv', row.names = FALSE)


### combinging outliers and removing
XCLUDE <- rbind(XCLUDE_UOE,XCLUDE_UEA)
# creating data-frame with control data removed
PMIfilt <- PMIdata[!(PMIdata$PPT %in% XCLUDE$PPT), ]

res_mediansF <- res_medians[!(res_medians$PPT %in% XCLUDE$PPT), ]
res_medians_allF <- res_medians_all[!(res_medians_all$PPT %in% XCLUDE$PPT), ]
res_meansF <- res_means[!(res_means$PPT %in% XCLUDE$PPT), ]
# saving filered data
setwd(anaPaths)
write.csv(PMIfilt, 'radialPMI-filtered.csv', row.names = FALSE)
write.csv(res_mediansF, 'radial-medians_filtered.csv', row.names = FALSE)
write.csv(res_meansF, 'radial-means_filtered.csv', row.names = FALSE)
write.csv(res_medians_allF, 'radial-medians-all_filtered.csv', row.names = FALSE)

### PLOT FILTERED PMI DATA
ggplot(PMIfilt, aes(x = SIDE, y = PMI, colour = SITE), position = position_dodge(.2)) + 
  geom_point(shape = 1, size = 4) +
  geom_line(aes(group = PPT), alpha = .5, size = .8) +
  stat_summary(aes(y = PMI, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1, size = 5, group = 1) +
  labs(title = 'Lateral Reaching', x = 'Side', y = 'Reaching error (mm)', 
       element_text(size = 14)) +
  facet_wrap(~DIAGNOSIS) +
  theme_bw() + theme(legend.position = 'right', text = element_text(size = 14),
                     strip.text.x = element_text(size = 12)) 

ggsave('radialPMI-filtered.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, width = 7, height = 4, path = anaPath)

## averaging PMI data across sides
meanFPMI <- summarySE(PMIfilt, measurevar = 'PMI', groupvar = c('DIAGNOSIS', 'DOM'),
                      na.rm = TRUE)
meanFPMI_all <- summarySE(PMIfilt, measurevar = 'PMI', groupvar = c('DIAGNOSIS'),
                          na.rm = TRUE)

#average across side
PMIfilt_av <- aggregate(PMI ~ PPT * SITE * DIAGNOSIS, mean, data = PMIfilt)
jitter <- position_jitter(width = 0.1, height = 0.1)

ggplot(PMIfilt_av, aes(x = DIAGNOSIS, y = PMI, colour = SITE)) + 
  geom_point(position = jitter, shape = 21, size = 3) +
  scale_colour_manual(values = c('grey40', 'black')) +
  stat_summary(aes(y = PMI, group = 1), fun.y = mean, colour = "black", 
               geom = 'point', shape = 3, stroke = 1, size = 4, group = 1) +
  labs(title = '', x = '', y = 'Reaching error (mm)', 
       element_text(size = 8)) +
  theme_bw() + theme(legend.position = 'none', text = element_text(size = 10),
                     strip.text.x = element_text(size = 8))

ggsave('radialPMI-filtered-av.png', plot = last_plot(), device = NULL, dpi = 300, 
       scale = 1, width = 3, height = 3, path = anaPath)

######## step 4: single case stats on filtered data ########
# create data-frames (controls and patients) with key information
td_summary <- summarySE(data = PMIfilt, measurevar = 'PMI', 
                        groupvars = c('DIAGNOSIS','DOM'), na.rm = TRUE)
# isolating control data for analysis
tdfilt_control <- td_summary[td_summary$DIAGNOSIS == 'HC', ]

#patient data for test of deficit
td_patient <- PMIdata[, c(1,3,4,11)] #extract from PMI because patient data not filtered
td_patient <- td_patient[td_patient$DIAGNOSIS != 'HC' ,] #removing controls
td_patient <- dcast(td_patient, PPT+DIAGNOSIS ~ DOM)
# NA values  = -1, so t-value is negative and can remove later
td_patient[is.na(td_patient$D), "D"] <- -1 

#patient data for test of deficit = same as above
# time for test of deficit! Calling on 'singcar' package developed by Jonathan Rittmo
# using Crawford's (1998) test of deficit
td_dom <- read.csv(text = 'PMI,TSTAT,PVALUE,PROP,ZCC,CI,LCI-T,HCI-T,LCI-P,HCI-P,DOM,PPT,DIAGNOSIS')
td_ndom <- read.csv(text = 'PMI,TSTAT,PVALUE,PROP,ZCC,CI,LCI-T,HCI-T,LCI-P,HCI-P,DOM,PPT,DIAGNOSIS')
for (l in 1:length(td_patient$PPT)){
  #left data first
  leftres <- TD(td_patient$ND[l], tdfilt_control$PMI[1], tdfilt_control$sd[1], 24, 
                alternative = 'greater', na.rm = FALSE)
  diff <- t(leftres$estimate)
  ltmp <- data.frame(td_patient$ND[l], leftres$statistic, leftres$p.value, 
                     diff[1], diff[2], t(leftres$interval), 'ND', check.names = FALSE) 
  ltmp$PPT <- td_patient$PPT[l]
  ltmp$DIAGNOSIS <- td_patient$DIAGNOSIS[l]
  td_ndom <- rbind(td_ndom, ltmp)
  #then right data
  rightres <- TD(td_patient$D[l], tdfilt_control$PMI[2], tdfilt_control$sd[2], 24, 
                 alternative = 'greater', na.rm = FALSE)
  diff <- t(rightres$estimate)
  rtmp <- data.frame(td_patient$D[l], rightres$statistic, rightres$p.value, 
                     diff[1], diff[2], t(rightres$interval), 'D', check.names = FALSE) 
  rtmp$PPT <- td_patient$PPT[l]
  rtmp$DIAGNOSIS <- td_patient$DIAGNOSIS[l]
  td_dom <- rbind(td_dom, rtmp)
}

# merging and renaming data-frames
tdfilt_results <- read.csv(text = 'PMI,TSTAT,PVALUE,ZCC,PROP,CI,LCI-T,HCI-T,LCI-P,HCI-P,DOM,PPT,DIAGNOSIS')
# changing names of td data-frames to match td-res
names(td_ndom) <- names(tdfilt_results)
names(td_dom) <- names(tdfilt_results)
tdfilt_results <- rbind(tdfilt_results, td_ndom, td_dom)
# remove original NA values (PMI = -1), where patients had no data
tdfilt_results <- tdfilt_results[tdfilt_results$PMI > -0.9 ,]
tdfilt_results$PPT <- factor(tdfilt_results$PPT)

# identifying patients with significant deficit
tdfilt_results$DEFICIT <- tdfilt_results$PVALUE < 0.025
# identifying patients with possible deficit
tdfilt_results$BL <- tdfilt_results$PVALUE < 0.05
# convert logicals into numerics, useful for binomial later
tdfilt_results$DEFICIT <- as.numeric(tdfilt_results$DEFICIT)
tdfilt_results$BL <- as.numeric(tdfilt_results$BL)

write.csv(tdfilt_results, 'radial-reaching_case-control.csv', row.names = FALSE)

## plotting p-values
ggplot(tdfilt_results, aes(x = DOM, y = PVALUE, group = PPT, colour = DIAGNOSIS)) +
  geom_point(size = 2, position = position_dodge(.2)) +
  geom_line(aes(group = PPT), size = 0.5, alpha = .5, position = position_dodge(.2)) + 
  scale_color_manual(values = c('black', 'grey50')) +
  geom_hline(yintercept = 0.025, color = 'red') +
  geom_hline(yintercept = 0.05, color = 'blue') +
  labs(title = 'Single case stats, filtered', x = 'Side', y = 'Probability of deficit', 
       element_text(size = 10)) +
  theme_bw()

ggsave('SCS_probability-scatter_filtered.png', plot = last_plot(), device = NULL, dpi = 300, 
       width = 6, height = 7, path = anaPath)

# Binomial statistics
# calculating the likelihood that inflated PMI occurs at above chance in each patient group
# with filtered data. Deficit and borderline deficit
pval <- .05
pval_Bl <- .1

## take any p < .025 for either side, re-run binomial on this data
# extracting p-value for each side (left and right) using td_results
td_sideF <- dcast(tdfilt_results, PPT+DIAGNOSIS ~ DOM, value.var = 'PVALUE')
td_sideF[is.na(td_sideF$D), "D"] <- 10
# finding defcits and borderline cases
td_sideF$DEFICIT <- (td_sideF$ND < 0.025 | td_sideF$D < 0.025)
td_sideF$BL <- (td_sideF$ND < 0.05 | td_sideF$D < 0.05)
# changing to numerical value
td_sideF$DEFICIT <- as.numeric(td_sideF$DEFICIT)
td_sideF$BL <- as.numeric(td_sideF$BL)

# splitting into MCI and AD
MCIPF <- td_sideF[td_sideF$DIAGNOSIS == 'MCI', ]
ADPF <- td_sideF[td_sideF$DIAGNOSIS == 'AD', ]

# binomial test of deficit
binMCIfilt <- binom.test(sum(MCIPF$DEFICIT), length(MCIPF$PPT), pval, alternative = 'greater')
print(binMCIfilt)
binADfilt <- binom.test(sum(ADPF$DEFICIT), length(ADPF$PPT), pval, alternative = 'greater')
print(binADfilt)
binALLfilt <- binom.test(sum(td_sideF$DEFICIT), length(td_sideF$PPT), pval, alternative = 'greater')
print(binALLfilt)

# binomial test of borderline deficit
binMCIfilt <- binom.test(sum(MCIPF$BL), length(MCIPF$PPT), pval_Bl, alternative = 'greater')
print(binMCIfilt)
binADfilt <- binom.test(sum(ADPF$BL), length(ADPF$PPT), pval_Bl, alternative = 'greater')
print(binADfilt)
binALLfilt <- binom.test(sum(td_sideF$BL), length(td_sideF$PPT), pval_Bl, alternative = 'greater')
print(binALLfilt)

## ANOVA ## 
# use PMIfilt data-frame to run between-subject ANOVA
# removing NA values for ANOVA - entire participant (not just side)
PMIanova <- PMIfilt[PMIfilt$PPT != 212, ]

# FULL ANOVA ON FILTERED DATA
FILT_ANOVA <- ezANOVA(
  data = PMIanova
  , dv = .(PMI)
  , wid = .(PPT)
  , within = .(DOM)
  , between = .(GRP)
  , type = 3
)

print(FILT_ANOVA)

### ANOVA BY ECCENTRICITY ###
# converting factor back to numeric keeping values
tmp <- as.numeric(levels(res_medians_allF$POSITION))[res_medians_allF$POSITION]

res_medians_allF$ECC <- abs(tmp)
ECCanova <- res_medians_allF[res_medians_allF$ECC != 400 ,]

ECCanova <- ECCanova[ECCanova$PPT != 212 ,]
ECCanova <- ECCanova[ECCanova$PPT != 101 ,]
ECCanova$ECC <- factor(ECCanova$ECC)

ECC_ANOVA <- ezANOVA(
  data = ECCanova
  , dv = .(AE)
  , wid = .(PPT)
  , within = .(DOM, VIEW, ECC)
  , between = .(GRP)
  , type = 3
)

print(ECC_ANOVA)
