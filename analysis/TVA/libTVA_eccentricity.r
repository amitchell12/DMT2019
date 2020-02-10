## Function that converts the OpenSesame combiTVA ouput into the stucture that is needed for parameter estimation in the libtva Matlab library

libTVA <- function(x, filename, filepath) {
  # Write trial number into the output file
  output.name <- paste0(filepath, filename, 'ALLtr', ".dat")

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
  # add an eccentricities column to letters
  letters$ECC <- tva_data$Eccentricity
  # code letter by by eccentricity
  for (l in 1:length(letters$ECC)){
    if (isTRUE(letters$ECC[l] == '1')){
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
   
  letters <- letters[, c(1:4,6:9)] 
  
  # create target and letter columns dependent on eccentricity
  ### reached here - keep going
  t1 <- t(letters)[t(positions == "top left")]
  t2 <- t(letters)[t(positions == "bot left")]
  t3 <- t(letters)[t(positions == "top right")]
  t4 <- t(letters)[t(positions == "bot right")]

  targets <- paste0(t1, t2, t3, t4)
    
  distractors <- rep("000000", length(tva_data$Block)) # No distractors in whole report condition
    
  
  ## Combine variables that are needed for the TVA parameter estimate
  whole_data <- cbind(
      condition = tva_data$Timing
      , letter_duration = tva_data$Letter_duration
      , targets
      , distractors
      , response = tva_data$response
      , ecc = tva_data$eEccentricity #adding eccentricity file
    )
    
    ## Write whole report data to file
    write.table(whole_data, file=output.name, quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t", append=TRUE)

   
}
 