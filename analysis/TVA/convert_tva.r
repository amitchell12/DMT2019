###########################################################################################################
# Before running this script, make sure the file "convert_tva.r" is in the same directory as this script  #
# and set the following paths correctly!                                                                  #
###########################################################################################################

#creating folder to sae data in
#dataPath <- ("S:/groups/DMT/data/") # Enter path to data
# Enter directory to save converted files to
#wPath <- "M:/GitHub/DMT2019/analysis/TVA/"
#setwd(dataPath)
#dir.create("formatted_TVA")
#outputPath <- ("S:/groups/DMT/data/formatted_TVA/")
#on mac
dataPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/data/'
wPath <- '/Users/alexandramitchell/Documents/git/DMT2019/analysis/TVA/'
setwd(dataPath)
dir.create('formatted_TVA')
outputPath <- '/Users/alexandramitchell/Documents/EDB_PostDoc/DMT2019/data/formatted_TVA/'
#back to script WD
setwd(wPath)


# Prepare files for MatLab fitting
if(readline("Have you removed all previously converted files? [y/n] ") == "n") {
  stop("Please delete or move the existing files.")
}

source("libTVA.r")
# listing all .csv files in data path
tva.files <- list.files(dataPath, ".csv", recursive = TRUE)
# finding those that just have 'TVA' in the title

setwd(dataPath)
for(i in tva.files) {
  if (isTRUE(substr(basename(i), 16, 16)=="w")){
    i.tva <- read.csv(paste0(dataPath, tva.files[which(tva.files == i)]))
    filename = gsub(".csv", "", i)
    libTVA(i.tva, basename(filename), filepath = outputPath)
  }
}

# Screen accuracy data
tva_accuracies <- data.frame(subject_nr = NA, accuracy = NA)
for (i in 1:length(tva.files)) {
  file = tva.files[i]
  if (isTRUE(substr(basename(file), 16, 16)=="w")){
    i_accuracies <- read.csv(paste0(dataPath, tva.files[i]))[c(10, 22)]
    tva_accuracies <- rbind(tva_accuracies, c(unique(i_accuracies$subject_nr), mean(as.numeric(as.character(i_accuracies$accumulated_accuracy)), na.rm = TRUE)))
  }
}

# removing first row (NA vals)
tva_accuracies <- tva_accuracies[-c(1), ]
boxplot(tva_accuracies$accuracy)
summary(tva_accuracies$accuracy)
# writing accuracy file
anaPath = 'S:/groups/DMT/analysis/TVA'
setwd(anaPath)
write.csv(tva_accuracies, 'TVAaccuracy.csv', row.names = FALSE)
