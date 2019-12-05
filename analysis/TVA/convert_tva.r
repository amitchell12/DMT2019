###########################################################################################################
# Before running this script, make sure the file "convert_tva.r" is in the same directory as this script  #
# and set the following paths correctly!                                                                  #
###########################################################################################################

#creating folder to sae data in
dataPath <- ("S:/groups/DMT/data/") # Enter path to data
wPath <- "M:/GitHub/DMT2019/analysis/TVA/"
# Enter directory to save converted files to
setwd(dataPath)
dir.create("formatted_TVA")
outputPath <- ("S:/groups/DMT/data/formatted_TVA/")
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

boxplot(tva_accuracies$accuracy)
summary(tva_accuracies$accuracy)
