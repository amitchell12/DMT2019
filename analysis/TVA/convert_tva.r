###########################################################################################################
# Before running this script, make sure the file "convert_tva.r" is in the same directory as this script  #
# and set the following paths correctly!                                                                  #
###########################################################################################################

#creating folder to sae data in
setwd("M:/Alex_Files/Experiments/DMT2019/DMT2019_rawdata/") # Enter path to 
dir.create("formatted_TVA")
#back to script WD
setwd("M:/Alex_Files/Experiments/DMT2019/DMT2019_rawdata/")

raw_data_path <- "M:/Alex_Files/Experiments/DMT2019/DMT2019_rawdata/" # Enter directory to raw data
output_path <- "M:/Alex_Files/Experiments/DMT2019/DMT2019_rawdata/formatted_TVA/" # Enter directory to save converted files to


# Prepare files for MatLab fitting
if(readline("Have you removed all previously converted files? [y/n] ") == "n") {
  stop("Please delete or move the existing files.")
}

source("libTVA.r")
tva.files <- list.files(raw_data_path, ".csv")

for(i in tva.files) {
  i.tva <- read.csv(paste0(raw_data_path, tva.files[which(tva.files == i)]))
  libTVA(i.tva, filename = gsub(".csv", "", i), filepath = output_path)
}

# Screen accuracy data
tva_accuracies <- data.frame(subject_nr = NA, accuracy = NA)
for (i in 1:length(tva.files)) {
  i_accuracies <- read.csv(paste0(raw_data_path, tva.files[i]))[seq(72, 468, 36), c("accumulated_accuracy", "subject_nr")]
  tva_accuracies <- rbind(tva_accuracies, c(unique(i_accuracies$subject_nr), mean(as.numeric(as.character(i_accuracies$accumulated_accuracy)))))
}

boxplot(tva.accuracies$accuracy)
summary(tva.accuracies$accuracy)