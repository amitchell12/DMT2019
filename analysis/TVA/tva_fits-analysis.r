##########################################################################################################
## ANALYSIS OF TVA FITS DATA MODELLED IN MATLAB ##
## AG.Mitchell 06.12.2019 ##

## first need to convert .txt files to readable .csv
#creating folder to sae data in
dataPath <- ("S:/groups/DMT/analysis/TVA/fits") # Enter path to data
anaPath <- "S:/groups/DMT/analysis/TVA/"
# Enter directory to save converted files to
setwd(dataPath)
dir.create("formatted_TVA")
outputPath <- ("S:/groups/DMT/data/formatted_TVA/")
#back to script WD
setwd(wPath)



