## Function that converts the OpenSesame combiTVA ouput into the stucture that is needed for parameter estimation in the libtva Matlab library

libTVA <- function(x, filename, filepath) {
  # Write trial number into the output file
  output.name <- paste0(filepath, filename, 'ALLtr', ".dat")
  output.name3 <- paste0(filepath, filename, 'ECC-3', '.dat')
  output.name9 <- paste0(filepath, filename, 'ECC-9', '.dat')
  
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
  positions <- tva_data[, paste0("posLet", c(1:4))]
  letters <- tva_data[, paste0("L", c(1:4))]
    
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
    )
    
    ## Write whole report data to file
    write.table(whole_data, file=output.name, quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t", append=TRUE)

    ##### small eccentricity set #####
    # by eccentricity
    tva_smallecc <- tva_data[tva_data$Eccentricity == 1 ,]
    exp_blocks <- subset(tva_smallecc, Block != 0) # Exclude practice trials
    write.table(length(exp_blocks$Block), file = output.name3, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "/t", append = FALSE)
    
    ## Sort targets by position (1-3 top left - bot left, 4-6 top right - bot right)
    positions <- tva_smallecc[, paste0("posLet", c(1:4))]
    letters <- tva_smallecc[, paste0("L", c(1:4))]
    
    t1 <- t(letters)[t(positions == "top left")]
    t2 <- t(letters)[t(positions == "bot left")]
    t3 <- t(letters)[t(positions == "top right")]
    t4 <- t(letters)[t(positions == "bot right")]
    
    targets <- paste0(t1, t2, t3, t4)
    
    distractors <- rep("000000", length(tva_smallecc$Block)) # No distractors in whole report condition
    
    
    
    # small ecc set
    whole_data_ecc3 <- cbind(
      condition = tva_smallecc$Timing
      , letter_duration = tva_smallecc$Letter_duration
      , targets
      , distractors
      , response = tva_smallecc$response
    )
    
    ## Write whole report data to file
    write.table(whole_data_ecc3, file=output.name3, quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t", append=TRUE)
    
    ##### large eccentricity set #####
    tva_largeecc <- tva_data[tva_data$Eccentricity == 2 ,]
    exp_blocks <- subset(tva_largeecc, Block != 0) # Exclude practice trials
    write.table(length(exp_blocks$Block), file = output.name9, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "/t", append = FALSE)
    
    
    ## Sort targets by position (1-3 top left - bot left, 4-6 top right - bot right)
    positions <- tva_largeecc[, paste0("posLet", c(1:4))]
    letters <- tva_largeecc[, paste0("L", c(1:4))]
    
    t1 <- t(letters)[t(positions == "top left")]
    t2 <- t(letters)[t(positions == "bot left")]
    t3 <- t(letters)[t(positions == "top right")]
    t4 <- t(letters)[t(positions == "bot right")]
    
    targets <- paste0(t1, t2, t3, t4)
    
    distractors <- rep("000000", length(tva_largeecc$Block)) # No distractors in whole report condition
    
    
    
    # small ecc set
    whole_data_ecc9 <- cbind(
      condition = tva_largeecc$Timing
      , letter_duration = tva_largeecc$Letter_duration
      , targets
      , distractors
      , response = tva_largeecc$response
    )
    
    ## Write whole report data to file
    write.table(whole_data_ecc9, file=output.name9, quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t", append=TRUE)
    
}
 