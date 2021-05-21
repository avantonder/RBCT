########################################################################################
### This is a R script containing analysis of all BASTA runs                         ###
### Helper functions have been taken or adapted from code used in Crispell 2019      ###
### https://github.com/JosephCrispell/GeneralTools/tree/master/WoodchesterPark/BASTA ###
########################################################################################

######################
### Load libraries ###
######################

library(ape)
library(tidyverse)
library(reshape2)
library(patchwork)

#################
### Functions ###
#################

##############################################
## Extract network names from log filenames ##
##############################################

getNetworkNames <- function(fileList){
  
  networkNames <- NULL
  
  for (i in 1:length(fileList)){
    networkNames[i] <- str_split(fileList, "_")[[i]][1]
  }
  
  networkNamesUnique <- unique(networkNames)
  
  return(networkNamesUnique)
}

#########################################
## Assign deme names to deme structure ##
#########################################

getDemeNamesForDemeStructure <- function(demeStructure, number=NULL){
  
  demeNames <- list(
    "2Deme"=c("Badger", "Bovine")
  )
  
  output = demeNames[[demeStructure]]
  if(is.null(number) == FALSE){
    output = demeNames[[demeStructure]][number + 1]
  }
  
  return(output)
}

####################################################
## If more than one replicate, combine log tables ##
####################################################

combineLogTables <- function(logTablesForReplicates){
  
  output <- logTablesForReplicates[[1]]
  if(length(logTablesForReplicates) > 1){
    for(i in 2:length(logTablesForReplicates)){
      output <- rbind(output, logTablesForReplicates[[i]])
    }
  }
  
  return(output)
}

#################################################################
## Convert backward migration rates to forward migration rates ##
#################################################################

calculateForwardMigrationRates <- function(logTable){
  
  # Get the names of the backward in time migration rate estimates
  migrationRateCols <- colnames(logTable)[
    grepl(colnames(logTable), pattern = "migModel.rateMatrix_")]
  
  # For each backward in time migration rate calculate the forward migration rate
  # FMR_ab = BMR_ba * (Nb / Na)
  #   MR: Migration rate (F - Forward, B - Backward)
  #   N: Effective population size
  #   Demes: a, b
  # Equation taken from second paragraph of Methods section in:
  # De Maio et al. 2015 - New routes to phylogeography ...
  #   BMR_ba = FMR_ab * (Na / Nb)
  #   ->  FMR_ab = BMR_ba / (Na / Nb)
  #     ->  FMR_ab = BMR_ba * (Nb / Na)
  #
  ### NOTE: Backward rates are multiplied by rate flag before being used ###
  # - Converts estimates to 0 when flag = 0 (i.e. rate turned off). If rate isn't likely, then they'll be a lot of zeros that will drag down estimate
  # - Set to NA when flag = 0
  
  for(backwardMigrationRateCol in migrationRateCols){
    
    # Get the demes involved
    parts <- strsplit(backwardMigrationRateCol, split="_")[[1]]
    a <- parts[2]
    b <- parts[3]
    
    #### Multiply the backward rates by the rate flag column ####
    backwardRate <- logTable[, backwardMigrationRateCol] * logTable[, paste("migModel.rateMatrixFlag_", a, "_", b, sep="")]
    
    # Get the estimate population sizes for a and b
    popASizes <- logTable[, paste("migModel.popSize_", a, sep="")]
    popBSizes <- logTable[, paste("migModel.popSize_", b, sep="")]
    
    # Calculate forward rate
    forwardMigrationRateCol <- paste("migModel.forwardRateMatrix_", b, "_", a, sep="")
    logTable[, forwardMigrationRateCol] <- backwardRate * (popASizes/popBSizes)
    
    # Convert the rates to NAs when flag set to zero
    logTable[, forwardMigrationRateCol][logTable[, paste("migModel.rateMatrixFlag_", a, "_", b, sep="")] == 0] <- NA
  }
  
  return(logTable)
}

#############################
## Read in BASTA log files ##
#############################

readInBASTALogTables <- function(networkNames, popEstimationTypes,
                                 clockEstimateTypes, path, burnInProp=0.1,
                                 nReplicates=NULL, ignoreIfFlagged=FALSE){
  
  # Store each of the log tables in a list
  logTables <- list()
  
  for(network in networkNames){
    
    for(popEstimationType in popEstimationTypes){
      
      for(clockEstimationType in clockEstimateTypes){
        
        # Build run defining prefix
        prefix <- paste(network, "_", "basta", "_", popEstimationType, "_",
                        clockEstimationType, sep="")
        
        # Initilialise a list to store the log tables for the current run's replicates
        logTablesForReplicates <- list()
        
        # Check if replicates available
        if(is.null(nReplicates) == FALSE){
          
          # Retrieve the data from each replicate
          for(replicate in 1:nReplicates){
            
            # Print progress information
            cat(paste("\rReading: ", prefix, "_", replicate, ".log\tReplicate: ", replicate,
                      "\t\t\t\t\t", sep=""))
            
            # Create file name
            file <- paste(path, prefix, "_", replicate, ".log", sep="")
            
            # Read in the file as table
            logTable <- read.table(file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
            
            # Replace "N" in table with NAs
            logTable[logTable == "N"] <- NA
            
            # Remove the burn-in
            burnIn <- round(burnInProp * nrow(logTable), digits=0)
            logTable <- logTable[burnIn:nrow(logTable), ]
            
            # Calculate the forward rates
            logTable <- calculateForwardMigrationRates(logTable)
            
            # Store the logTable
            logTablesForReplicates[[as.character(replicate)]] <- logTable
          }
          
          # Store the tables as a single log table
          logTables[[prefix]] <- combineLogTables(logTablesForReplicates)
          
        }else{
          
          # Print progress information
          cat(paste("\rReading: ", prefix, ".log\t\t\t\t\t\t", sep=""))
          
          # Create file name
          file <- paste(path, prefix, ".log", sep="")
          
          # Read in the file as table
          logTable <- read.table(file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
          
          # Replace "N" in table with NAs
          logTable[logTable == "N"] <- NA
          
          # Remove the burn-in
          burnIn <- round(burnInProp * nrow(logTable), digits=0)
          logTable <- logTable[burnIn:nrow(logTable), ]
          
          # Calculate the forward rates
          logTable <- calculateForwardMigrationRates(logTable)
          
          # Store the tables as a single log table
          logTables[[prefix]] <- logTable
        }
      }
    }
  }
  cat("\rFinished reading in log tables...\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\n")
  
  return(logTables)
}

########################################
## Calculate ESS values for each plot ##
########################################

calculateEffectiveSampleSize <- function(posteriorSample){
  
  # Calculate the mean of the sample
  sampleMean <- mean(posteriorSample)
  nSamples <- length(posteriorSample)
  
  # Create an array to store the correlation values between a[index] vs. a[index + lag]
  autoCorrelationValues <- c()
  
  for(lag in 1:(nSamples - 1)){
    
    # Calculate the auto-correlation value for the current lag period
    # Taken from: http://www.itl.nist.gov/div898/handbook/eda/section3/autocopl.htm
    a <- posteriorSample[1:(nSamples - lag)] - sampleMean
    b <- posteriorSample[(1:(nSamples - lag) + lag)] - sampleMean
    
    Ch <- (1/nSamples) * sum(a * b)
    C0 <- sum((posteriorSample - sampleMean)^2) / nSamples
    autoCorrelationValues[lag] <- Ch / C0
    
    # Check if auto-correlation has dropped below zero
    if(lag != 1 && 
       (autoCorrelationValues[lag - 1] > 0 && autoCorrelationValues[lag] <= 0 || 
        autoCorrelationValues[lag - 1] < 0 && autoCorrelationValues[lag] >= 0)){
      break;
    }
    
    # Monitor progress
    #if(lag %% 1000 == 0){
    #  cat(paste("Calculate correlation based upon gap of ", lag, ". Max gap = ", (nSamples - 1), "\n", sep=""))
    #}
  }
  
  # Calculate the Effective Sample Size
  # Taken from: http://people.duke.edu/~ccc14/sta-663-2016/16C_PyMC3.html
  ess <- nSamples / (1 + 2 * (sum(autoCorrelationValues)))
  
  return(ess);
}

########################################
## Plot ESS values for each parameter ##
########################################

plotParameterESSValues <- function(logTable, colNamesToPlot){
  
  essValues <- rep(NA, length(colNamesToPlot))
  
  for(i in 1:length(colNamesToPlot)){
    
    # Get the posterior values from the current column
    values <- as.numeric(logTable[, colNamesToPlot[i]])
    
    # Remove transition rate values when rate flag is set to zero
    if(grepl(colNamesToPlot[i], pattern="migModel.rateMatrix") == TRUE){
      
      # Get the deme numbers
      parts = strsplit(colNamesToPlot[i], split="_")[[1]]
      a <- parts[2]
      b <- parts[3]
      
      # Build the rate Flag column name
      flagCol <- paste("migModel.rateMatrixFlag_", a, "_", b, sep="")
      
      # Remove values where rate flag == 0
      values <- values[logTable[, flagCol] != 0]
    }
    
    # Remove NAs, if present
    values <- values[is.na(values) == FALSE]
    
    # Note and skip those where the value is always the same
    range <- range(values)
    if(range[1] == range[2]){
      #cat(paste("Parameter \"", colNamesToPlot[i],
      #          "\" always has single value: ", range[1], "\n", sep=""))
      next
    }
    
    if(length(essValues) > 1){
      essValues[i] <- calculateEffectiveSampleSize(values)
    }else{
      essValues[i] <- length(essValues[i])
    }
  }
  
  #Create dataframe for plotting
  
  essValuesMelt <- melt(essValues)
  
  extractedESSvalues <- cbind(colNamesToPlot, essValuesMelt)
  
  # Create plot
  
  extractedESSvaluesPlot <- ggplot(data = extractedESSvalues, aes(x = colNamesToPlot,y = value, fill = colNamesToPlot)) + 
    geom_bar(stat = "identity") +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
    ylab("Effective Sample Size") +
    xlab("Posterior parameter") +
    scale_fill_discrete(name = "Posterior parameter") +
    ggtitle("Parameter Effective Sample Sizes") +
    geom_hline(yintercept = 100, colour = "red", linetype = "dashed") +
    theme_bw()
  
  print(extractedESSvaluesPlot)
  
}

################################################################
## Plot posterior support for each deme as root (2 Deme only) ##
################################################################

plotPosteriorSupportForEachDemeAsRoot <- function(logTable, demeStructure){
  
  RootPosteriorSupport <- logTable$treePrior.rootColor
  
  RootPosteriorSupport[RootPosteriorSupport == 0] <- "Badger"
  
  RootPosteriorSupport[RootPosteriorSupport == 1] <- "Cow"
  
  RootPosteriorSupportMelt <- melt(RootPosteriorSupport)
  
  RootPosteriorSupportPlot <- ggplot(data = RootPosteriorSupportMelt, aes(as.factor(value), fill = as.factor(value))) + 
    geom_bar() +
    xlab("Host") + ylab("Count") + theme_bw() +
    scale_fill_discrete(name = "Host") +
    ggtitle("Assignment of Demes to Root State")
  
  print(RootPosteriorSupportPlot)
  
}

##################################################
## Plot population size histogram for each deme ##
##################################################

plotPopulationSizes <- function(logTable, demeStructure){
  
  # Get the population size distributions
  
  popSizeEstimates <- list()
  
  for(col in colnames(logTable)){
    
    if(grepl(col, pattern="popSize") == FALSE){
      next
    }
    
    popSizeEstimates[[strsplit(col, split="_")[[1]][2]]] <- logTable[, col]
  }
  
  # Get the deme names
  demeNames <- getDemeNamesForDemeStructure(demeStructure)
  
  names(popSizeEstimates) <- demeNames
  
  popSizeEstimatesCombined <- bind_rows(popSizeEstimates)
  
  popSizeEstimatesCombinedMelt <- melt(popSizeEstimatesCombined)
  
  popSizeEstimatePlot <- ggplot(data= popSizeEstimatesCombinedMelt, aes(x = value, group = variable, fill = variable)) + 
    geom_density(adjust=1.5, alpha = 0.5) + xlab("Effective size") + ylab("Density") + theme_bw() +
    scale_fill_discrete(name = "Host") +
    ggtitle("Deme Effective Population Sizes")
  
  print(popSizeEstimatePlot)  
}

####################################################
## Calculate and plot substitution rate estimates ##
####################################################

examineSubstitutionRateEstimates <- function(logTable, genomeSize){
  
  # Get the substitution rate estimates note that stored differently between strict and relaxed
  rateEstimates <- c()
  
  if(length(which(grepl(colnames(logTable), pattern="mutationRate") == TRUE)) == 1){
    
    rateEstimates <- logTable$mutationRate
  }else{
    
    rateEstimates <- logTable$ucedMean
  }
  
  Likelihood <- as.numeric(cut(logTable$posterior, breaks = 10))
  
  substitutionRatePlot <- ggplot() + geom_point(aes(x = rateEstimates * genomeSize, y = logTable$tree.height, 
                                                    colour = Likelihood)) +
    theme_bw() + 
    xlab("Substitution Rate (per Genome per Year)") + 
    ylab("Root Height (years)") + 
    ggtitle("Substitution Rate versus Root Height")
  
  print(substitutionRatePlot)
  
  return(rateEstimates * genomeSize)
}

###########################################
## Summarise substitution rate estimates ##
###########################################

summariseDistribution <- function(distribution){
  
  # Calculate the median
  median <- median(distribution, na.rm=TRUE)
  
  # Calculate the lower and upper bounds
  quantiles <- quantile(distribution, probs=c(0.025, 0.975), na.rm=TRUE)
  
  cat(paste0("median = ", median, " (lower = ", quantiles[1], ", upper = ", quantiles[2], ")"))
}

#################################################################
## Calculate migration rates between demes for arrow thickness ##
#################################################################

getArrowRates <- function(logTable){
  
  # Initialise a list to store the arrow weights
  arrowRates <- list()
  
  # Get the column names
  colNames <- colnames(logTable)
  
  # Examine each column
  for(col in colNames){
    
    # Ignore all columns except the forward rate columns
    if(grepl(x=col, pattern="forward") == FALSE){
      next
    }
    
    # Skip rate if never estimate e.d. cattle-outer -> badger-inner
    if(length(unique(logTable[, col])) <= 2){
      next
    }
    
    # Get the directional information (i.e. 0_1, 2_0)
    parts <- strsplit(col, split="_")[[1]]
    direction <- paste(parts[length(parts) - 1], "_", parts[length(parts)], sep="")
    
    # Store a summary statistic for the rate distribution
    arrowRates[[direction]] <- logTable[, col]
  }
  
  return(arrowRates)
}

################################
## Extract values from a list ##
################################

getValues <- function(list){
  values <- c()
  keys <- names(list)
  for(i in 1:length(keys)){
    values[i] <- list[[keys[i]]]
  }
  
  return(values)
}

######################
## Normalise values ##
######################

divideValuesInListByMax <- function(arrowRates){
  max <- max(getValues(arrowRates), na.rm=TRUE)
  output <- timesValuesInListByValue(arrowRates, 1/max)
  
  return(output)
}

###################################################
## Multiply migration rates by normalised values ##
###################################################

timesValuesInListByValue <- function(list, value){
  
  output <- list()
  for(key in names(list)){
    output[[key]] <- list[[key]] * value
  }
  
  return(output)
}

##########################
## Create an empty plot ##
##########################

createEmptyPlot <- function(){
  par(mar=c(0,0,0,0))
  plot(x=NULL, y=NULL, xlim=c(0,1), ylim=c(0,1), xaxt="n", yaxt="n", bty="n",
       ylab="", xlab="")
}

###############################################################
## Plot migration rates between cattle and badgers as arrows ##
###############################################################

plotMigrationRates <- function(logTable, demeStructure, code, arrowFactor){
  
  # Get the migration rates
  output <- getArrowRates(logTable)
  
  # Calculate the median for each migration rate being estimate
  medians <- list()
  for(rate in names(output)[names(output) != "AICM"]){
    medians[[rate]] <- median(output[[rate]], na.rm=TRUE)
  }
  
  # Normalise those rates to vary between 0 and MAX (arrow factor)
  migrationRates <- divideValuesInListByMax(medians)
  migrationRates <- timesValuesInListByValue(migrationRates, arrowFactor)
  
  # Check for "NaN" or NA values
  for(key in names(migrationRates)){
    if(is.nan(migrationRates[[key]]) == TRUE || is.na(migrationRates[[key]]) == TRUE){
      migrationRates[[key]] <- 0
    }
  }
  
  demeStructure <- "2Deme"
  
  # Create empty plot
  createEmptyPlot()
  
  # Get deme names and assign colours
  demeNames <- getDemeNamesForDemeStructure(demeStructure)
  demeColours <- rep("black", length(demeNames))
  demeColours[grepl(demeNames, pattern="badger")] <- "red"
  demeColours[grepl(demeNames, pattern="cow")] <- "blue"
  
  # Add labels
  x <- c(0.1, 0.9)
  y <- c(0.5, 0.5)
  text(x=x, y=y, 
       labels=demeNames,
       col=demeColours)
  
  # badger -> cow
  if(migrationRates[["0_1"]] != 0){
    arrows(x0=x[1]+0.15, x1=x[2]-0.15, y0=y[1]-0.15, y1=y[2]-0.15,
           code=code, lwd=migrationRates[["0_1"]])
  }
  # cow -> badger
  if(migrationRates[["1_0"]] != 0){
    arrows(x0=x[2]-0.15, x1=x[1]+0.15, y0=y[2]+0.15, y1=y[1]+0.15,
           code=code, lwd=migrationRates[["1_0"]])
  }
  
  return(output)
}

#####################################
## Plot migration rate  posteriors ##
#####################################

plotMigrationRatePosteriors <- function(logTable, demeStructure){
  
  # Get the deme names
  demeNames <- getDemeNamesForDemeStructure(demeStructure)
  
  # Get the forward migration rate estimates and their associated flags
  rateEstimates <- list()
  rateFlags <- list()
  for(col in colnames(logTable)){
    
    if(grepl(col, pattern="migModel.forwardRateMatrix_") == TRUE){
      
      # Get the demes involved
      parts <- strsplit(col, split="_")[[1]]
      a <- parts[2]
      b <- parts[3]
      
      # Skip rate if never estimate e.d. cattle-outer -> badger-inner
      if(length(unique(logTable[, col])) == 1){
        next
      }
      
      # Build migration rate name
      rateName <- paste(demeNames[as.numeric(a) + 1], "->", demeNames[as.numeric(b) + 1])
      
      # Store the posterior distribution
      rateEstimates[[rateName]] <- logTable[, col]
      
      # Find the rate flag values for this rate (note that flags will be for backward rate from b to a)
      rateFlags[[rateName]] <- logTable[, paste("migModel.rateMatrixFlag_", b, "_", a, sep="")]
    }
  }
  
  # Plot a histogram for each migration rate posterior
  for(key in names(rateEstimates)){
    
    values <- rateEstimates[[key]]
    values[is.na(values)] <- 0
    
    valuesDF <- as.data.frame(values)
    
    migrationRatePlot <- ggplot(valuesDF) + geom_histogram(aes(values)) + 
      theme_bw() + 
      ylab("Count") +
      xlab("Migration rate posterior") +
      ggtitle(key)
    
    print(migrationRatePlot)
    
  }
  
  migrationRatesMeans <- NULL
  
  for(i in 1:length(rateFlags)){
    
    migrationRatesMeans[[i]] <- mean(rateFlags[[i]])
    
  }  
  
  names(migrationRatesMeans) <- names(rateFlags)
  
  migrationRatesMeansMelt <- melt(migrationRatesMeans)
  
  migrationRatesMeansPlot <- ggplot(data = migrationRatesMeansMelt) + geom_point(aes(x = row.names(migrationRatesMeansMelt), y = value)) +
    xlab("") +
    ylab("Proportion MCMC steps flag = ON") +
    ggtitle("Forward migration rate posterior flags") +
    theme_bw()
  
  print(migrationRatesMeansPlot)
  
}

#########################################
## Plot MCMC traces for each parameter ##
#########################################

plotParameterTraces <- function(logTable, colNamesToPlot){
  
  # Define a grid of plots
  nPlots <- ceiling(sqrt(length(colNamesToPlot)))
  
  # Set the plotting window dimensions
  par(mfrow=c(nPlots, nPlots))
  
  # Set the margins
  par(mar=c(0, 0, 2, 0))
  
  # Plot a trace for each parameter
  for(col in colNamesToPlot){
    plot(logTable[, col], type="l", xlab="", xaxt="n", ylab="", yaxt="n", bty="n",
         main=col, cex.main=0.5)
  }
  
  # Reset the margins
  par(mar=c(5.1, 4.1, 4.1, 2.1))
  
  # Reset the plotting window dimensions
  par(mfrow=c(1,1))
}

###########################
## Calculate AICM values ##
###########################

calculateAICM <- function(logLikelihoodValues, nBootstraps){
  # Calculation taken from:
  # Baele et al. 2012 - Improving the accuracy of demographic and molecular clock
  # model comparison while accomodating phylogenetic uncertainty
  # Equation 10: AICM = k - 2m
  #   k = effective number of parameters = 2 * variance of the posterior log likehoods
  #   m = mean of the posterior log likelihoods
  aicm <- (2 * var(logLikelihoodValues)) - (2 * mean(logLikelihoodValues))
  
  # Conduct some bootstrapping
  bootstrapAICMs <- c()
  for(i in 1:nBootstraps){
    
    sample <- sample(logLikelihoodValues, size=length(logLikelihoodValues), replace=TRUE)
    bootstrapAICMs[i] <- (2 * var(sample)) - (2 * mean(sample))
  }
  
  # Summarise the bootstraps
  quantiles <- quantile(bootstrapAICMs, probs=c(0.025, 0.975))
  output <- c(aicm, quantiles[[1]], quantiles[[2]])
  
  return(output)
}

getValues <- function(list){
  values <- c()
  keys <- names(list)
  for(i in 1:length(keys)){
    values[i] <- list[[keys[i]]]
  }
  
  return(values)
}

####################################################################################
## Produce summaries for log files for each network and model. Output to PDF file ##
####################################################################################

summarisePosteriorLogTables <- function(path, logTables, code=2, arrowFactor, nBootstraps,
                                        genomeSize){
  
  # Get analysis names
  analyses <- names(logTables)
  
  # Initialise a list to store the migration rate estimates and AICM
  migrationRateEstimates <- list()
  
  # Create directory for summary plots
  dir.create(paste(path, "SummaryPlots", "/", sep=""), showWarnings = FALSE)
  
  # Examine each analysis
  for(analysis in analyses){
    
    # Open a pdf file
    file <- paste(path, "SummaryPlots", "/", analysis, "_ResultsSummary.pdf",
                  sep="")
    pdf(file)
    
    # Get the deme structure (only using 2 deme model)
    demeStructure <- "2Deme"
    
    # Note progress
    cat(paste("\rExamining log table for: ", analysis, "\t\t\t\t\t\t\t", sep=""))
    
    # Get the log table
    logTable <- logTables[[analysis]]
    
    # Note which parameters to examine posterior support for
    colsToCalculateESS <- colnames(logTable)[
      grepl(x=colnames(logTable), pattern="Sample|rateMatrixFlag|forward|Count") == FALSE]
    
    # Plot the ESS values of each parameter estimated
    plotParameterESSValues(logTable, colsToCalculateESS)
    
    # Plot the posterior support for each deme as source
    plotPosteriorSupportForEachDemeAsRoot(logTable, demeStructure)
    
    # Plot the population size estimates
    plotPopulationSizes(logTable, demeStructure)
    
    # Examine the substitution rate estimates
    rateEstimates <- examineSubstitutionRateEstimates(logTable, genomeSize)
    
    if(grepl(analysis, pattern = "equal")){
      cat("\n", analysis, "\n")
      summariseDistribution(rateEstimates)
      cat("\n")
    }
    
    # Produce a migration rate estimation figure - weight by rate flags
    # Diagrams designed with code = 2 (FORWARDS) in mind
    migrationRateEstimates[[analysis]] <- 
      plotMigrationRates(logTable, demeStructure, code, arrowFactor)
    
    # Plot the forward migration rate posterior distributions and their flags
    plotMigrationRatePosteriors(logTable, demeStructure)
    
    # Plot the parameter traces
    plotParameterTraces(logTable, colsToCalculateESS)
    
    # Calculate the aicm
    migrationRateEstimates[[analysis]][["AICM"]] <- 
      calculateAICM(logTable$treeLikelihood1, nBootstraps)
    
    
    # Close the pdf file
    dev.off()
  }
  cat("\rFinished examining log tables...\t\t\t\t\t\t\n")
  
  # Reset the margins
  
  return(migrationRateEstimates)
}

##########################################
## Calculate AICM scores for each model ##
##########################################

getAICMScores <- function(migrationRateEstimates){
  
  # Initialise a list to store the analysis AICM scores
  scores <- list()
  
  # Examine each analysis
  for(key in names(migrationRateEstimates)){
    
    scores[[key]] <- migrationRateEstimates[[key]]$AICM
  }
  
  return(scores)
}

############################
## Plot AICM model scores ##
############################

plotModelAICMScores <- function(migrationRateEstimates, nBootstraps){
  
  # Get the deme structure names
  analyses <- names(migrationRateEstimates)
  
  # Initialise an array to store the shortened analyses names
  network <- c()
  names <- c()
  
  # Get the AICM values
  aicmScores <- c()
  bootstrapMins <- c()
  bootstrapMaxs <- c()
  
  for(i in 1:length(analyses)){
    
    values <- migrationRateEstimates[[analyses[i]]][["AICM"]]
    
    aicmScores[i] <- values[1]
    bootstrapMins[i] <- values[2]
    bootstrapMaxs[i] <- values[3]
    
    parts <- strsplit(analyses[i], split="_")[[1]]
    network[i] <- parts[1]
    names[i] <- paste(parts[3], parts[4], sep=" ")
  }
  
  AICMdf <- as.data.frame(cbind(network, names, aicmScores, bootstrapMins, bootstrapMaxs), stringsAsFactors = FALSE)
  
  # Colour minimum AICM score red
  
  AICMdfColour <- AICMdf %>%
    group_by(network) %>%
    mutate(colour = (min(aicmScores) == aicmScores))
  
  AICMplot <- ggplot(data = AICMdfColour) + 
    geom_point(aes(x = factor(names), y = as.numeric(aicmScores), color = colour)) + 
    geom_errorbar(aes(x = factor(names), ymin = as.numeric(bootstrapMins), ymax = as.numeric(bootstrapMaxs), color = colour), 
                  width = 0.2, position = position_dodge(0.9)) + 
    theme_bw() + ylab("AICM score") + xlab("BASTA model") + ylim(min(as.numeric(bootstrapMins)) - 1, 
                                                                 max(as.numeric(bootstrapMaxs)) + 1) + 
    theme(axis.text.y = element_text(size = 15), axis.text.x = element_text(size = 15, angle = 90), 
          axis.title.y = element_text(size = 15), axis.title.x = element_text(size = 15)) +
    theme(legend.position = "none") +
    facet_wrap(~network, ncol = 4) +
    scale_color_manual(values = c("black", "red"))
  
  return(AICMplot)
  
}

####################################
## Read in transition count files ##
####################################

readTransitionCountFiles <- function(migrationRateEstimates, path, nReplicates=NULL){
  
  transitionCounts <- list()
  
  analyses <- names(migrationRateEstimates)
  
  for (analysis in analyses){
    
    transitionCountsForReplicates <- list()
    
    prefix <- paste(analysis, "_", "TransitionCounts", sep="")
    
    for(replicate in 1:nReplicates){
      
      # Print progress information
      cat(paste("\rReading: ", prefix, "_", replicate, ".txt\tReplicate: ", replicate,
                "\t\t\t\t\t", sep=""))
      
      # Create file name
      file <- paste(path, prefix, "_", replicate, ".txt", sep="")
      
      # Read in the file as table
      transitionCount <- read.table(file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
      
      transitionCount$Replicate <- replicate
      
      transitionCount <- transitionCount[,c(ncol(transitionCount),1:ncol(transitionCount)-1)]
      
      # Store the logTable
      transitionCountsForReplicates[[as.character(replicate)]] <- transitionCount
    }
    
    transitionCounts[[analysis]] <- combineLogTables(transitionCountsForReplicates)
  }
  
  return(transitionCounts)
}

##########################################
## Remove burnin from transition counts ##
##########################################

removeBurnIn <- function(transitionCounts, burnInProp=0.1, nReplicates=3){
  
  analyses <- names(transitionCounts)
  
  # Examine each analysis
  for(analysis in analyses){
    
    # Examine each replicate
    for(rep in 1:nReplicates){
      
      # Get the indices assocated with the current analysis
      rows <- which(transitionCounts[[analysis]]$Replicate == rep)
      
      # Note the order of the rows - based on the posterior sample
      order <- order(transitionCounts[[analysis]][rows, "Sample"], decreasing=FALSE)
      
      # Note the first burn-in indices
      burnIn <- round(burnInProp * length(rows), digits=0)
      
      indicesToRemove <- rows[order][1:burnIn]
      
      # Remove the burn-in rows from the transition counts
      transitionCounts[[analysis]] <- transitionCounts[[analysis]][-indicesToRemove, ]
    }
  }
  
  return(transitionCounts)
}

###########################################################################
## Plot estimated transition rates for each network using the best model ##
###########################################################################

calculateMeanEstimatedTransitionRatesBetweenCattleAndBadgerPopulationsWeightedByAICM <-
  function(migrationRateEstimates){
    
    # Get the analyses names
    analyses <- names(migrationRateEstimates)
    names <- c()
    
    # Initialise vectors to store the median and range of the total rates between cattle and badgers
    modelBadgerToCowRateMedians <- c()
    modelBadgerToCowRateLower <- c()
    modelBadgerToCowRateUpper <- c()
    modelCowToBadgerRateMedians <- c()
    modelCowToBadgerRateLower <- c()
    modelCowToBadgerRateUpper <- c()
    
    # Initialise vectors to summarise the distribution ratio: badger-to-cattle / cattle-to-cattle
    modelInterSpeciesRatioMedians <- c()
    modelInterSpeciesRatioLower <- c()
    modelInterSpeciesRatioUpper <- c()
    
    # Initialise a vectors to store the mean of the rate flags
    modelBadgerToCowFlagMeans <- c()
    modelCowToBadgerFlagMeans <- c()
    
    # Examine each of the different model structures
    for(i in 1:length(analyses)){
      
      # Initialise arrays to store the rates from badgers to cattle and vice versa
      sumRatesBadgerToCow <- rep(0, length(migrationRateEstimates[[analyses[i]]][[1]]))
      sumRatesCowToBadger <- rep(0, length(migrationRateEstimates[[analyses[i]]][[1]]))
      
      # Initialise arrays to store the rate flags from badgers to cattle and vice versa
      sumFlagsBadgerToCow <- rep(0, length(migrationRateEstimates[[analyses[i]]][[1]]))
      sumFlagsCowToBadger <- rep(0, length(migrationRateEstimates[[analyses[i]]][[1]]))
      
      # Initialise counts for the number of badger-to-cow and cow-to-badger rates for the current model
      nBadgerToCow <- 0
      nCowToBadger <- 0
      
      # Get the demeStructure
      parts <- strsplit(analyses[i], "_")[[1]]
      demeStructure <- "2Deme"
      names[i] <- parts[1]
      
      # Examine each rate
      for(key in names(migrationRateEstimates[[analyses[i]]])){
        
        # Ignore AICM and zero rates (rates never estimate e.g. outer-cattle -> inner-badger)
        if(key == "AICM"){
          next
        }
        
        # Get the migration rate values
        values <- migrationRateEstimates[[analyses[i]]][[key]]
        
        # Create the flags (1 when values not NA, 0 otherwise)
        flags <- ifelse(is.na(values), 0, 1)
        
        # Replace the NA values with zero - so they aren't carried into the sums
        values[is.na(values)] <- 0
        
        # Split the key into its deme numbers
        demeNumbers <- as.numeric(strsplit(key, split="_")[[1]])
        
        # Check if current rate is between badger and cattle populations
        if(grepl(getDemeNamesForDemeStructure(demeStructure, demeNumbers[1]),
                 pattern="Badger") == TRUE &&
           grepl(getDemeNamesForDemeStructure(demeStructure, demeNumbers[2]),
                 pattern="Bovine") == TRUE){
          
          # Increment the count
          nBadgerToCow <- nBadgerToCow + 1
          
          # Add the flags for the current rate to the growing sum
          sumFlagsBadgerToCow <- sumFlagsBadgerToCow + flags
          
          # Add the current rate estimates to the growing sum
          sumRatesBadgerToCow <- sumRatesBadgerToCow + values
          
        }else if(grepl(getDemeNamesForDemeStructure(demeStructure, demeNumbers[1]),
                       pattern="Bovine") == TRUE &&
                 grepl(getDemeNamesForDemeStructure(demeStructure, demeNumbers[2]),
                       pattern="Badger") == TRUE){
          
          # Increment the count
          nCowToBadger <- nCowToBadger + 1
          
          # Add the flags for the current rate to the growing sum
          sumFlagsCowToBadger <- sumFlagsCowToBadger + flags
          
          # Add the current rate estimates to the growing sum
          sumRatesCowToBadger <- sumRatesCowToBadger + values
        }
      }
      
      # Calculate the ratio of the inter-species rate sums
      interSpeciesTransitionRateRatio <- sumRatesBadgerToCow / sumRatesCowToBadger
      interSpeciesTransitionRateRatio <- interSpeciesTransitionRateRatio[interSpeciesTransitionRateRatio != Inf]
      modelInterSpeciesRatioMedians[i] <- median(interSpeciesTransitionRateRatio)
      quantiles <- quantile(interSpeciesTransitionRateRatio, probs=c(0.025, 0.975))
      modelInterSpeciesRatioLower[i] <- quantiles[[1]]
      modelInterSpeciesRatioUpper[i] <- quantiles[[2]]
      
      # Remove any values that are exactly zero
      # - These will result when flag=0 across badger-to-cattle/cattle-to-badger rates estimate
      sumRatesBadgerToCow <- sumRatesBadgerToCow[sumRatesBadgerToCow != 0]
      sumRatesCowToBadger <- sumRatesCowToBadger[sumRatesCowToBadger != 0]
      
      # Store summary statistics of these posterior sums
      modelBadgerToCowRateMedians[i] <- median(sumRatesBadgerToCow)
      quantiles <- quantile(sumRatesBadgerToCow, probs=c(0.025, 0.975))
      modelBadgerToCowRateLower[i] <- quantiles[[1]]
      modelBadgerToCowRateUpper[i] <- quantiles[[2]]
      
      modelCowToBadgerRateMedians[i] <- median(sumRatesCowToBadger)
      quantiles <- quantile(sumRatesCowToBadger, probs=c(0.025, 0.975))
      modelCowToBadgerRateLower[i] <- quantiles[[1]]
      modelCowToBadgerRateUpper[i] <- quantiles[[2]]
      
      # Divide the sums of the flags by the number of rates
      flagsBadgerToCow <- sumFlagsBadgerToCow / nBadgerToCow
      flagsCowToBadger <- sumFlagsCowToBadger / nCowToBadger
      
      # Store summary statistics of these posterior sums
      modelBadgerToCowFlagMeans[i] <- mean(flagsBadgerToCow)
      modelCowToBadgerFlagMeans[i] <- mean(flagsCowToBadger)
    }
    
    ## Create dataframe of all calculated values
    
    transmissionRatesDF <- as.data.frame(cbind(names, modelBadgerToCowRateLower, modelBadgerToCowRateMedians, 
                                               modelBadgerToCowRateUpper, modelCowToBadgerRateLower, modelCowToBadgerRateMedians, 
                                               modelCowToBadgerRateUpper, modelBadgerToCowFlagMeans,
                                               modelCowToBadgerFlagMeans), stringsAsFactors = FALSE)
    
    #return(transmissionRatesDF)
    
    EstimatedTransitionRatesBetweenCattleAndBadgerPlot <- ggplot(data = transmissionRatesDF) +
      geom_point(aes(x = names, y = as.numeric(modelBadgerToCowRateMedians)), color = "darkred", 
                 position = position_nudge(x = -0.1)) + 
      geom_errorbar(aes(x = names, ymin = as.numeric(modelBadgerToCowRateLower), 
                        ymax = as.numeric(modelBadgerToCowRateUpper)), color = "darkred", 
                    width = 0.2, position = position_nudge(x = -0.1)) + 
      geom_point(aes(x = names, y = as.numeric(modelCowToBadgerRateMedians)), color = "darkblue",
                 position = position_nudge(x = 0.1)) + 
      geom_errorbar(aes(x = names, ymin = as.numeric(modelCowToBadgerRateLower), 
                        ymax = as.numeric(modelCowToBadgerRateUpper)), color = "darkblue", 
                    width = 0.2, position = position_nudge(x = 0.1)) +
      geom_text(aes(x = names, y = as.numeric(modelBadgerToCowRateUpper), 
                    label = round(as.numeric(modelBadgerToCowFlagMeans), digits = 2)), color = "darkred",
                position = position_nudge(x = -0.1, y = 0.1), size = 3) +
      geom_text(aes(x = names, y = as.numeric(modelCowToBadgerRateUpper), 
                    label = round(as.numeric(modelCowToBadgerFlagMeans), digits = 2)), color = "darkblue",
                position = position_nudge(x = 0.1, y = 0.1), size = 3) +
      theme_bw() + 
      ylab("Per lineage transition rate per year") + 
      xlab("Transmission cluster")
    
    return(EstimatedTransitionRatesBetweenCattleAndBadgerPlot)
  }

plotSummaryOfTransitionCountsBasedOnPosteriorTrees <- function(transitionCountEstimates){
  
  analyses <- names(transitionCountEstimates)
  
  # Note the columns of interest
  columns <- c("Count_Badger.Badger", "Count_Bovine.Bovine", "Count_Badger.Bovine", "Count_Bovine.Badger")
  
  Median <- NULL
  FirstQuantile <- NULL
  ThirdQuantile <- NULL
  Combined <- NULL
  allTransitions <- NULL
  
  
  for (analysis in analyses){
    
    for(column in columns){
      
      # Calculate the quantiles for the current rate
      
      Median[[column]] <- median(transitionCountEstimates[[analysis]][, column])
      FirstQuantile[[column]] <- quantile(transitionCountEstimates[[analysis]][, column], 0.025)
      ThirdQuantile[[column]] <- quantile(transitionCountEstimates[[analysis]][, column], 0.975)
      Combined <- rbind(Median, FirstQuantile, ThirdQuantile)
      Combined <- as.data.frame(t(Combined))
      Combined$Transition <- row.names(Combined)
      
      allTransitions[[analysis]] <- Combined
      
    }
    
  }
  
  allTransitionsDF <- bind_rows(allTransitions, .id = "model")
  
  allTransitionsDF <- allTransitionsDF %>%
    mutate(model2 = model) %>% 
    separate(model2, c("network", NA, NA, NA), sep="_")  %>% 
    mutate(Transition2 = str_replace(Transition, "Count_", "")) %>% 
    mutate(Transition3 = str_replace(Transition2, "[.]", "-"))
  
  allTransitionsPlot <- ggplot(data = transform(allTransitionsDF, 
                                                Transition4 = factor(Transition3, 
                                                                     levels = c("Badger-Badger", "Bovine-Bovine",                                                     
                                                                                "Badger-Bovine", "Bovine-Badger")),
                                                network = factor(network, 
                                                                 levels = c(1,2,3,4,5,6,7,8,9,10,11,12)))) +
    geom_point(aes(x = Transition4, y = Median, color = Transition4), size = 3) +
    geom_errorbar(aes(x = Transition4, ymin = FirstQuantile, 
                      ymax = ThirdQuantile, color = Transition4), width = 0.3) +
    theme_bw() +
    ylab("Number of transitions") +
    xlab("") +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(color = "Host transition") + 
    theme(axis.title=element_text(size=15)) + 
    theme(axis.text.x = element_text(size = 15)) + 
    theme(axis.text.y = element_text(size = 15)) +
    theme(legend.position = "none") +
    facet_wrap(~network, ncol = 4) + 
    theme(strip.text.x = element_text(size = 15, color = "black"),
          strip.background = element_rect(color="black", fill="white", size=1, linetype="solid"))
  
  
  return(allTransitionsPlot)
}

#####################################
## Set location of BASTA log files ##
#####################################

BASTAlogPath <- "~/All_BASTA/"

##########################
## List BASTA log files ##
##########################

BASTAlogFiles <- list.files(pattern = "*.log")

###################################
## Extract network names as list ##
###################################

BASTAlogNetworkNames <- getNetworkNames(BASTAlogFiles)

################################
## Note deme structure to use ##
################################

demeStructures <- list(
  "2Deme" = 2)

#############################################
## List population size estimation options ##
#############################################

popEstimationTypes <- c("varying", "equal")

##################################
## List the clock model options ##
##################################

clockEstimateTypes <- c("relaxed", "strict") 

#############################
## Read in BASTA log files ##
#############################

logTables <- readInBASTALogTables(BASTAlogNetworkNames, popEstimationTypes,
                                  clockEstimateTypes, BASTAlogPath, nReplicates=3)

#######################################################
## Specify number of bootstraps for calculating AICM ##
#######################################################

nBootstraps <- 1000

############################################################
## Define the size of the genome; Here masked genome size ##
############################################################

genomeSize <- 4135167

######################################
## Summarise results from log files ##
######################################

migrationRateEstimates <- summarisePosteriorLogTables(BASTAlogPath, logTables, code=2,
                                                      arrowFactor=20, nBootstraps,
                                                      genomeSize)

##################################
## Note each model's AICM score ##
##################################

aicmScores <- getAICMScores(migrationRateEstimates)

###################################
## Examine the model likelihoods ##
###################################

aicmOrderPlot <- plotModelAICMScores(migrationRateEstimates, nBootstraps)

################################################
## Collapse list of AICM score into single df ##
################################################

aicmScoresDFAll <- bind_rows(aicmScores, .id = "model")

aicmScoresDFAll <- melt(aicmScoresDFAll[1,])

aicmScoresDFAll <- aicmScoresDFAll %>% 
  select(variable, value) %>%
  mutate(model = variable) %>% 
  separate(model, c("network", NA, NA, NA), sep="_")


#############################################
## Use equal strict model for all networks ##
#############################################

bestAICMscores <- aicmScoresDFAll %>%
  filter(grepl('equal_strict', variable))

bestAICMscoresList <- bestAICMscores[["variable"]]

#########################################################
## Get the migration rate estimates for the best model ##
#########################################################

bestModelmigrationRateEstimates <- NULL

for (i in seq_along(bestAICMscoresList)){
  
  name <- bestAICMscoresList[i]
  
  bestModelmigrationRateEstimates[[name]] <- migrationRateEstimates[[name]]
}

bestModelmigrationRateEstimatesCompact <- compact(bestModelmigrationRateEstimates)

names(bestModelmigrationRateEstimatesCompact) <- bestAICMscoresList

#######################################################################
## Plot transition rates between cattle and badgers for all clusters ##
#######################################################################

TransitionRatesBetweenCattleAndBadgerPlot <- calculateMeanEstimatedTransitionRatesBetweenCattleAndBadgerPopulationsWeightedByAICM(bestModelmigrationRateEstimatesCompact)

##########################################
## Calculate aggregate transition rates ##
##########################################

BadgerToBovineTransitions <- NULL

BovineToBadgerTransitions <- NULL

for (i in 1:length(bestModelmigrationRateEstimatesCompact)){
  
  BadgerToBovineTransitions[[i]] <- bestModelmigrationRateEstimatesCompact[[i]][2]
  BovineToBadgerTransitions[[i]] <- bestModelmigrationRateEstimatesCompact[[i]][1]
}

names(BadgerToBovineTransitions) <- BASTAlogNetworkNames

BadgerToBovineTransitionsDF <- bind_rows(BadgerToBovineTransitions, .id = "network")

names(BovineToBadgerTransitions) <- BASTAlogNetworkNames

BovineToBadgerTransitionsDF <- bind_rows(BovineToBadgerTransitions, .id = "network")

CombinedTransitionsDF <- as_tibble(cbind(BadgerToBovineTransitionsDF, BovineToBadgerTransitionsDF[,2]))

CombinedTransitionsDF <- CombinedTransitionsDF %>% 
  mutate(flags_0_1 = ifelse(is.na(`0_1`), 0, 1)) %>% 
  mutate(flags_1_0 = ifelse(is.na(`1_0`), 0, 1))

CombinedTransitionsDF[is.na(CombinedTransitionsDF)] <- 0

CombinedTransitionsDFfiltered <- CombinedTransitionsDF %>%
  mutate(interSpeciesTransitionRateRatio = `0_1`/`1_0`) %>% 
  filter(interSpeciesTransitionRateRatio != Inf) %>% 
  filter(`0_1` != 0) %>% 
  filter(`1_0` != 0)

modelBadgerToCowRateMedians <- median(CombinedTransitionsDFfiltered$`0_1`)
modelCowToBadgerRateMedians <- median(CombinedTransitionsDFfiltered$`1_0`)
modelBadgerToCowRateQuantiles <- quantile(CombinedTransitionsDFfiltered$`0_1`, probs=c(0.025, 0.975))
modelCowToBadgerRateQuantiles <- quantile(CombinedTransitionsDFfiltered$`1_0`, probs=c(0.025, 0.975))
modelBadgerToCowRateLower <- modelBadgerToCowRateQuantiles[[1]]
modelBadgerToCowRateUpper <- modelBadgerToCowRateQuantiles[[2]]
modelCowToBadgerRateLower <- modelCowToBadgerRateQuantiles[[1]]
modelCowToBadgerRateUpper <- modelCowToBadgerRateQuantiles[[2]]
modelBadgerToCowFlagMeans <- mean(CombinedTransitionsDFfiltered$flags_0_1)
modelCowToBadgerFlagMeans <- mean(CombinedTransitionsDFfiltered$flags_1_0)

names <- "All"

AllTransmissionRatesDF <- as.data.frame(cbind(names, modelBadgerToCowRateLower, modelBadgerToCowRateMedians, 
                                              modelBadgerToCowRateUpper, modelCowToBadgerRateLower, modelCowToBadgerRateMedians, 
                                              modelCowToBadgerRateUpper, modelBadgerToCowFlagMeans,
                                              modelCowToBadgerFlagMeans), stringsAsFactors = FALSE)

##############################################################################################
## Combine transmissionRatesDF from TransitionRatesBetweenCattleAndBadgerPlot function with ##
## AllTransmissionRatesDF containing aggregate results for all transition rates             ##
##############################################################################################

CombinedTransmissionRatesDF <- rbind(transmissionRatesDF, AllTransmissionRatesDF)

#######################################################################################
## Plot transition rates for all clusters and aggregate transition rates (Figure 2E) ##
#######################################################################################

EstimatedTransitionRatesBetweenCattleAndBadgerPlotAggregate <- ggplot(data = transform(CombinedTransmissionRatesDF,
                                                                                       names = factor(names, levels = c("All",1,2,3,4,5,6,7,8,9,10,11,12)))) +
  geom_point(aes(x = names, y = as.numeric(modelBadgerToCowRateMedians)), color = "#ff0000", 
             position = position_nudge(x = -0.1)) + 
  geom_errorbar(aes(x = names, ymin = as.numeric(modelBadgerToCowRateLower), 
                    ymax = as.numeric(modelBadgerToCowRateUpper)), color = "#ff0000", 
                width = 0.2, position = position_nudge(x = -0.1)) + 
  geom_point(aes(x = names, y = as.numeric(modelCowToBadgerRateMedians)), color = "#0000cc",
             position = position_nudge(x = 0.1)) + 
  geom_errorbar(aes(x = names, ymin = as.numeric(modelCowToBadgerRateLower), 
                    ymax = as.numeric(modelCowToBadgerRateUpper)), color = "#0000cc", 
                width = 0.2, position = position_nudge(x = 0.1)) +
  geom_text(aes(x = names, y = as.numeric(modelBadgerToCowRateUpper), 
                label = round(as.numeric(modelBadgerToCowFlagMeans), digits = 2)), color = "#ff0000",
            position = position_nudge(x = -0.1, y = 0.005), size = 4) +
  geom_text(aes(x = names, y = as.numeric(modelCowToBadgerRateUpper), 
                label = round(as.numeric(modelCowToBadgerFlagMeans), digits = 2)), color = "#0000cc",
            position = position_nudge(x = 0.1, y = 0.005), size = 4) +
  theme_bw() + 
  ylab("Per lineage transition rate per year") + 
  xlab("Transmission cluster") + 
  theme(axis.title=element_text(size=15)) + 
  theme(axis.text.x = element_text(size = 15)) + 
  theme(axis.text.y = element_text(size = 15))

#####################################################################################################
### Examine the transition counts recorded on the posterior tree distributions from each analysis ###
#####################################################################################################

################################################################################
## Read in the transition counts for all the trees                            ##
## Counts created by Java: MyWork:examineBASTAPosterior:CountTransitions.java ##
################################################################################

transitionCounts <- readTransitionCountFiles(migrationRateEstimates, BASTAlogPath, nReplicates=3)

##################################################
## Remove the burn-in period from each analysis ##
##################################################

transitionCountsWithoutBurnIn <- removeBurnIn(transitionCounts, burnInProp=0.1, nReplicates=3)

##################################################
## Get the transition counts for the best model ##
##################################################

bestModelTransitionCounts <- NULL

for (i in seq_along(bestAICMscoresList)){
  
  name <- bestAICMscoresList[i]
  
  bestModelTransitionCounts[[name]] <- transitionCountsWithoutBurnIn[[name]]
}

bestModelTransitionCountsCompact <- compact(bestModelTransitionCounts)

names(bestModelTransitionCountsCompact) <- bestAICMscoresList

#########################
## Read in host counts ##
#########################

hostCounts <- read_csv("network_animals.csv")

bestAICMscoresHostCounts <- bestAICMscores %>% 
  left_join(hostCounts, by = c("network" = "network"))

#############################################################
## Plot transition counts from posterior trees (Figure 2F) ##
#############################################################

transitionsCountsPlot <- plotSummaryOfTransitionCountsBasedOnPosteriorTrees(bestModelTransitionCountsCompact)

