#############################
## Load required libraries ##
#############################

library(ape)
library(reshape2)
library(TransPhylo)
library(coda)
library(parallel)
library(tidyverse)

###################################
## Define getInfectTime function ##
###################################

getInfectTime <- function(record, k){
  mytimes = vapply(1:length(record), function(x) {
        tt = extractTTree(record[[x]]$ctree)
        ii = which(tt$nam == k)
        return(tt$ttree[ii, 1])
    }, FUN.VALUE = 1)
  return(mytimes)
}

###########################################################################
## Read in BEAST tree (rooted/unrooted). The tree labels need to be in   ##
## the format: sample_name-decimal_date e.g. sample1-2009.345            ##
###########################################################################

network1_original_tree1 <- ggtree::read.beast("network1_masked_snps_dated_strict_constant_184_MCC.tree")

##############
## Set seed ##
##############

set.seed(1)

#########################################################################
## Prepare dtr tree by removing any negative branch lengths to prevent ##
## later Transphylo error                                              ##
#########################################################################

############################
## Remove multifurcations ##
############################

network1_dtr_no_multifurcations1 <- ape::multi2di(network1_original_tree1@phylo)

network1_dtr_no_multifurcations1$edge.length <- pmax(network1_dtr_no_multifurcations1$edge.length, 1/365)

#################################################
## Convert BEAST tree to TransPhylo input tree ##
#################################################

network1_ptree1 <- TransPhylo::ptreeFromPhylo(network1_dtr_no_multifurcations1, dateLastSample = 2006.926027)

####################################
## Define TransPhylo run function ##
####################################

runTransphylo <- function(i){
  TransPhylo::inferTTree(network1_ptree1, mcmcIterations = 1000000, w.shape = 1.6, w.scale = 3.5, ws.shape = 1.6, ws.scale = 3.5, thinning = 50, dateT = 2006.927)
}

##################################
## Run TransPhylo in triplicate ##
##################################

network1_Transphylo_runs1 <- parallel::mclapply(1:3, runTransphylo)

####################################################################
## Convert Transphylo MCMC runs to coda format                    ## 
## @ burnin (proportion of MCMC output to be discarded as burnin) ##
####################################################################

network1_Transphylo_coda_results1 <- parallel::mclapply(network1_Transphylo_runs1, function(x) TransPhylo::convertToCoda(x, burnin=0.1))

##############################
## Add Coda outputs to list ##
##############################

network1_Transphylo_coda_list1 <- coda::mcmc.list(network1_Transphylo_coda_results1[[1]], network1_Transphylo_coda_results1[[2]], network1_Transphylo_coda_results1[[3]])

#########################################################################################
## Assess whether Transphylo runs have converged by running Gelman And Rubin's         ##
## Convergence Diagnostic on all three coda objects.  Ideally results should be < 1.05 ##
#########################################################################################

network1_Transphylo_converg_test1 <- coda::gelman.diag(network1_Transphylo_coda_list1, confidence = 0.95, transform = FALSE, autoburnin = FALSE, multivariate = TRUE)

###################################################################################
## Assuming MCMC runs have converged, use first MCMC run for subsequent analyses ##
###################################################################################

#################################################################################################
## Calculate effective population size values (ESS) for first MCMC run. Values should be > 100 ##
#################################################################################################

network1_Transphylo_ESS1 <- coda::effectiveSize(network1_Transphylo_coda_results1[[1]])

#######################################
## Save TransPhylo run to Rdata file ##
#######################################

save.image(file = "network1_transphylo_beast_1.Rdata")

###########################################
## Calculate consensus transmission tree ##
###########################################

network1_ttree1 <- TransPhylo::consTTree(network1_Transphylo_runs1[[1]])

network1_tt1 <- extractTTree(network1_Transphylo_runs1[[1]][[1]]$ctree)

###################
## Remove burnin ##
###################

network1_Transphylo_run_post_burnin1 <- network1_Transphylo_runs1[[1]][max(1, round(length(network1_Transphylo_runs1[[1]]) * 0.5)):length(network1_Transphylo_runs1[[1]])]

###########################################
## Extract sample names by removing date ##
###########################################

network1_samples1 <- network1_tt1$nam

network1_names1 <- NULL

for (i in 1:length(network1_samples1)){
  network1_names1[i] <- paste(str_split(network1_samples1, "_")[[i]][1], str_split(network1_samples1, "_")[[i]][2],sep = '_')
}

###########################################
## Extract sample dates by removing name ##
###########################################

network1Samplingdates1 <- NULL

for (i in 1:length(network1_samples1)){
  network1Samplingdates1[i] <- strsplit(network1_samples1, "_")[[i]][3]
}

#################################################
## Create data frame of sample names and dates ##
#################################################

network1Samplingdf1 <- as.data.frame(cbind(network1_names1,network1Samplingdates1), stringsAsFactors = FALSE)

############################################################
## Extract estimated times between infection and sampling ##
############################################################

network1_infectionTimes1 <- lapply(network1_samples1, function(x) {
  getInfectTime(network1_Transphylo_run_post_burnin1, x)
})

network1_infectionTimesdf1 <- as.data.frame((do.call("cbind", lapply(network1_infectionTimes1, ts))))

colnames(network1_infectionTimesdf1) <- network1_names1

######################################
## Save infection times to csv file ##
######################################

write_csv(network1_infectionTimesdf1, "network1_beast_infection_times_1.csv")

########################
## Save to Rdata file ##
########################

save.image(file = "network1_transphylo_beast_1.Rdata")


