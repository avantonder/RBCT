library(TransPhylo)
library(tidyverse)
library(igraph)
library(tidygraph)

######################
## Read in metadata ##
######################

network1Metadata <- read_csv("network1_metadata.csv")

################################
## Load TransPhylo Rdata file ##
################################

network1TransPhyloResults <- load("network1_transphylo_beast_1.Rdata")

################################################################
## Extract ESS values for TransPhylo run and save to csv file ##
################################################################

write_csv(as.data.frame(t(network1_Transphylo_ESS1)), "network1_Transphylo_ESS.csv")

##################################################################################
## Extract sampled and unsampled cases from TransPhylo run and save to csv file ##
##################################################################################

network1incidentCases <- getIncidentCases(network1_Transphylo_runs1[[1]])

network1incidentCasesdf <- data.frame(Time=c(network1incidentCases$Time,network1incidentCases$Time),Cases=c(network1incidentCases$sampledCases,network1incidentCases$unsampCases), IsSampled=c(rep("Yes",length(network1incidentCases$sampledCases)),rep("No",length(network1incidentCases$unsampCases))))

write_csv(network1incidentCasesdf, "network1_incident_cases.csv")

##########################################
## Calculate Transmission probabilities ##
##########################################

network1TransmissionProb <- computeMatWIW(network1_Transphylo_runs1[[1]])

#########################################################
## Rename transmission probability matrix column names ##
#########################################################

network1TransmissionProbColNames <- NULL

for (i in 1:length(colnames(network1TransmissionProb))){
  network1TransmissionProbColNames[i] <- paste(str_split(colnames(network1TransmissionProb), "_")[[i]][1], str_split(colnames(network1TransmissionProb), "_")[[i]][2],sep = '_')
}

colnames(network1TransmissionProb) <- network1TransmissionProbColNames

#########################################################
## Rename transmission probability matrix row names ##
#########################################################

network1TransmissionProbRowNames <- NULL

for (i in 1:length(row.names(network1TransmissionProb))){
  network1TransmissionProbRowNames[i] <- paste(str_split(row.names(network1TransmissionProb), "_")[[i]][1], str_split(row.names(network1TransmissionProb), "_")[[i]][2],sep = '_')
}

row.names(network1TransmissionProb) <- network1TransmissionProbRowNames

##################################################################################################
## Convert transmission probability matrix to new df showing transmission probabilities between ##
## each pair of isolates                                                                        ##
##################################################################################################

network1TransmissionProbDecon <- data.frame(t(combn(row.names(network1TransmissionProb),2)), prob = t(network1TransmissionProb)[lower.tri(network1TransmissionProb)])
colnames(network1TransmissionProbDecon) <- c("Taxon1", "Taxon2", "prob")

#####################################################
## Save out transmission probabilities as csv file ##
#####################################################

write_csv(network1TransmissionProbDecon, "network1_transmission_probabilities.csv")

###############################
## Extract transmission tree ##
###############################

network1TransmissionTree <- extractTTree(network1_Transphylo_runs1[[1]][[1]]$ctree)

network1TransmissionInfo <- network1TransmissionTree$ttree

network1TransmissionSimpleWiw <- cbind(network1TransmissionInfo[,3],1:length(network1TransmissionInfo[,1]))

####################################
## Remove dates from sample names ##
####################################

network1TransmissionNames <- NULL

for (i in 1:length(network1TransmissionTree$nam)){
  network1TransmissionNames[i] <- paste(str_split(network1TransmissionTree$nam, "_")[[i]][1], str_split(network1TransmissionTree$nam, "_")[[i]][2],sep = '_')
}

##############################
## Create nodes for network ##
##############################

network1TransmissionNodes <- data.frame(id = 1:nrow(network1TransmissionInfo), label=c(network1TransmissionNames,(length(network1TransmissionTree$nam)+1) : nrow(network1TransmissionInfo)), stringsAsFactors = FALSE)

######################################
## Add host and trial id to node df ##
######################################

network1TransmissionNodes$host <- network1Metadata$Host[match(network1TransmissionNodes$label, network1Metadata$lane_id)]
network1TransmissionNodes$trial_id <- network1Metadata$trial_id[match(network1TransmissionNodes$label, network1Metadata$lane_id)]

##########################################
## Add unsampled '0' isolate to node df ##
##########################################

network1TransmissionNodes[is.na(network1TransmissionNodes)] <- "Unsampled"

network1TransmissionNodesnrows <- nrow(network1TransmissionNodes) + 1

network1TransmissionNodes[network1TransmissionNodesnrows,]$id <- 0
network1TransmissionNodes[network1TransmissionNodesnrows,]$label <- "0"
network1TransmissionNodes[network1TransmissionNodesnrows,]$host <- "Unsampled"
network1TransmissionNodes[network1TransmissionNodesnrows,]$trial_id <- "Unsampled"

########################################################
## Save out nodes of transmission network as csv file ##
########################################################

write_csv(network1TransmissionNodes, "network1TransmissionNodes.csv")

##############################
## Create edges for network ##
##############################

network1TransmissionEdges <- data.frame(from = network1TransmissionSimpleWiw[,1], to = network1TransmissionSimpleWiw[,2], arrows="to")

########################################################
## Save out egdes of transmission network as csv file ##
########################################################

write_csv(network1TransmissionEdges, "network1TransmissionEdges.csv")