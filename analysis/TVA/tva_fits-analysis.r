##########################################################################################################
## ANALYSIS OF TVA FITS DATA MODELLED IN MATLAB ##
## AG.Mitchell 06.12.2019 ##

##### file info #####
#creating folder to sae data in
dataPath <- ("S:/groups/DMT/analysis/TVA/fits/") # Enter path to data
anaPath <- "S:/groups/DMT/analysis/TVA/"
# Enter directory to save converted files to
setwd(dataPath)

# list all files in working directory
txt_filelist <- list.files(dataPath, ".txt")

for (i in 1:length(txt_filelist)) {
  FILE = read.table(file=txt_filelist[i], header = TRUE, sep = '\t')
  
  write.table(
    FILE,
    file=paste0(dataPath,
                sub(".txt","",txt_filelist[i]),".csv"),
    row.names=FALSE,
    quote = FALSE,
    sep=",")
}
 
