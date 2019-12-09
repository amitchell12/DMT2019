##########################################################################################################
## ANALYSIS OF TVA FITS DATA MODELLED IN MATLAB ##
## AG.Mitchell 06.12.2019 ##

##### file info #####
#creating folder to sae data in
dataPath <- 'S:/groups/DMT/data/formatted_TVA/'
fitPath <- ("S:/groups/DMT/analysis/TVA/fits/") # Enter path to data
anaPath <- "S:/groups/DMT/analysis/TVA/"
# Enter directory to save converted files to
setwd(fitPath)

# list all files in working directory
txt_filelist <- list.files(fitPath, ".txt")

# converting to .csv files
for (i in 1:length(txt_filelist)) {
  FILE = read.table(file=txt_filelist[i], header = TRUE, sep = '\t')
  
  write.table(
    FILE,
    file=paste0(fitPath,
                sub(".txt","",txt_filelist[i]),".csv"),
    row.names=FALSE,
    quote = FALSE,
    sep=",")
}
 
##### extracting key data #####
# first fing t0 for each participant
### can't get this to work - try and fix
dat_filelist <- list.files(dataPath, ".dat")
for (i in dat_filelist){
  if (substr(basename(dat_filelist[i]), -4, -4)){
    
  }
}