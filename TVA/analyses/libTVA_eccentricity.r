## Function that converts the OpenSesame combiTVA ouput into the stucture that is needed for parameter estimation in the libtva Matlab library

libTVA <- function(x, filename, filepath) {
  # Write trial number into the output file
  output.name <- paste0(filepath, filename, ".dat")

  exp_blocks <- subset(x, Block != 0) # Exclude practice trials
  write.table(length(exp_blocks$Block), file = output.name, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "/t", append = FALSE)
  
  # Whole report data
  tva_data <- subset(exp_blocks)
  tva_data$Letter_duration <- tva_data$Letter_duration # Correct to represent actual presentation times
  tva_data[, paste0("L", 1:4)] <- apply(tva_data[, paste0("L", 1:4)], 2, as.character)
    
  ## Reformat responses (this should be fixed in OpenSesame!)
  tva_data$response <- toupper(tva_data$response) # Convert to upper case
  tva_data$response <- gsub(" ", "", tva_data$response)
  tva_data$response[tva_data$response == "NONE" | tva_data$response == ""] <- "-" # Code missing as "-"
    
  ## Sort targets by position (1-3 top left - bot left, 4-6 top right - bot right)
  ###### unsure about how to work this- keep going
  positions <- tva_data[, paste0("posLet", c(1:4))]
  letters <- tva_data[, paste0("L", c(1:4))]
  positions[] <- lapply(positions, as.character) #changing from factor to character
  # add an eccentricities column to letters and positions
  letters$ECC <- tva_data$Eccentricity
  positions$ECC <- tva_data$Eccentricity
  # code letter by by eccentricity
  for (l in 1:length(letters$ECC)){
    if (isTRUE(letters$ECC[l] == '2')){
      letters$L6[l] = letters$L1[l]
      letters$L7[l] = letters$L2[l]
      letters$L8[l] = letters$L3[l]
      letters$L9[l] = letters$L4[l]
      letters[l, c(1:4)] = 0
    }
    else {
      letters$L6[l] = 0
      letters$L7[l] = 0
      letters$L8[l] = 0
      letters$L9[l] = 0
      }
      
  }
  
  
  ##### this doesn't work - figure out why
  # do the same for positions
  for (l in 1:length(positions$ECC)){
    if (isTRUE(positions$ECC[l] == '2')){
      positions$posLet6[l] = positions$posLet1[l]
      positions$posLet7[l] = positions$posLet2[l] 
      positions$posLet8[l] = positions$posLet3[l] 
      positions$posLet9[l] = positions$posLet4[l] 
      positions[l, c(1:4)] = 0
    }
    else {
      positions$posLet6[l] = 0
      positions$posLet7[l] = 0
      positions$posLet8[l] = 0
      positions$posLet9[l] = 0
    }
    
  }
   
  letters <- letters[, c(1:4,6:9)] 
  positions <- positions[, c(1:4,6:9)]
  
  # create target and letter columns dependent on eccentricity
  # add target columns 5:8
t1 <- t(letters)[t(positions == "top left")]
t2 <- t(letters)[t(positions == "bot left")]
t3 <- t(letters)[t(positions == "top right")]
t4 <- t(letters)[t(positions == "bot right")]
  
for (l in 1:length(tva_data$Eccentricity)){
  if (isTRUE(tva_data$Eccentricity[l] == 1)){
    targets[l] <- paste0(t1[l], t2[l], t3[l], t4[l], 0, 0, 0, 0)
  }
  else {
    targets[l] <- paste0(0, 0, 0, 0, t1[l], t2[l], t3[l], t4[l])
  }
}


    
distractors <- rep("00000000", length(tva_data$Block)) # No distractors in whole report condition
    
  
  ## Combine variables that are needed for the TVA parameter estimate
  whole_data <- cbind(
      condition = tva_data$Timing
      , letter_duration = tva_data$Letter_duration
      , targets
      , distractors
      , response = tva_data$response
    )
    
    ## Write whole report data to file
    write.table(whole_data, file=output.name, quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t", append=TRUE)

   
}
 