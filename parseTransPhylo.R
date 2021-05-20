###############################
### Load required libraries ###
###############################

library(tidyverse)
library(reshape2)
library(igraph)
library(tidygraph)
library(ggraph)

#############################
### Set working directory ###
#############################

allTransPhyloDir <- "~/All_Transphylo/"
setwd(allTransPhyloDir)

######################
## Read in metadata ##
######################

allMetadata <- read_csv("all_metadata_110520.csv")

##################################
## Read in sensible cluster IDs ##
##################################

newNetworkIDs <- read_csv("new_network_ids.csv")

#################################
## Read in pairsnp output file ##
#################################

M_bovisPWsnps <- read_csv("ERADbTB_no_Bennett_no_modern_new_151119_masked_snps.csv")

###########################################
## Convert pairsnp output into dataframe ##
###########################################

M_bovisPWsnpsDecon <- data.frame( t(combn(names(M_bovisPWsnps),2)), dist=t(M_bovisPWsnps)[lower.tri(M_bovisPWsnps)] )
colnames(M_bovisPWsnpsDecon) <- c("Taxon1", "Taxon2", "dist")
M_bovisPWsnpsDecon$Taxon1 <- gsub('# ', "",M_bovisPWsnpsDecon$Taxon1)

########################
## List ESS csv files ##
########################

ESSFiles <- list.files(pattern = "*_ESS.csv")

###################################
## Extract cluster names as list ##
###################################

ESSFilesNames <- gsub('_Transphylo_ESS.csv', '', ESSFiles)

#######################
## Read in csv files ##
#######################

ESSFilesRead <- lapply(ESSFiles, function(x){
  read_csv(x)
})

names(ESSFilesRead) <- ESSFilesNames

###############################################
## Collapse list of ESS score into single df ##
###############################################

ESSFilesDFAll <- bind_rows(ESSFilesRead, .id = "network_id")

ESSFilesDFAll <- ESSFilesDFAll %>% 
  select(-off.p)

######################################
## Plot ESS scores for each cluster ##
######################################

ESSFilesDFAllMelt <- melt(ESSFilesDFAll)

ESSPlot <- ggplot(data = ESSFilesDFAllMelt, aes(x = variable,y = value, fill = variable)) + 
  geom_bar(stat = "identity") +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
  ylab("Effective Sample Size") +
  xlab("Posterior parameter") +
  scale_fill_discrete(name = "Posterior parameter") +
  theme_bw() +
  coord_flip() +
  facet_wrap(~network_id, ncol = 4) +
  geom_hline(yintercept = 100, colour = "red", linetype = "dashed")

###################################
## List infection time csv files ##
###################################

infectFiles <- list.files(pattern = "*_1.csv")

###################################
## Extract cluster names as list ##
###################################

infectFilesNames <- gsub('_beast_infection_times_1.csv', '', infectFiles)

#######################
## Read in csv files ##
#######################

infectFilesRead <- lapply(infectFiles, function(x){
  read_csv(x)
})

names(infectFilesRead) <- infectFilesNames

#####################################################
## Calculate median infection time for each sample ##
#####################################################

InfectTimesMedian <- lapply(infectFilesRead, function(x){
  sapply(x, median)
})

###################################################################  
## Create dataframe with sample names and median infection times ##
###################################################################

InfectTimesMedianDF <- lapply(InfectTimesMedian, function(x){
  as_tibble(cbind(names(x), x))
})

############################################################
## Add cluster id and host to each infection times tibble ##
############################################################

for (i in  1:length(InfectTimesMedianDF)){
  names(InfectTimesMedianDF[[i]]) <- c("sample", "median")
  #InfectTimesMedianDF[[i]]$network_id <- infectFilesNames[i]
  InfectTimesMedianDF[[i]]$host <- allMetadata$host[match(InfectTimesMedianDF[[i]]$sample, allMetadata$lane_id)]
  InfectTimesMedianDF[[i]]$sampling_date <- allMetadata$decimal_date[match(InfectTimesMedianDF[[i]]$sample, allMetadata$lane_id)]
  InfectTimesMedianDF[[i]]$infect_time <- InfectTimesMedianDF[[i]]$sampling_date - as.numeric(InfectTimesMedianDF[[i]]$median)
}

#####################################################
## Collapse list of infection times into single df ##
#####################################################

InfectTimesMedianDFAll <- bind_rows(InfectTimesMedianDF, .id = "network_id")

InfectTimesMedianDFAll <- InfectTimesMedianDFAll %>% 
  left_join(newNetworkIDs, by = "network_id")

####################################################
## Calculate median infection times for each host ##
####################################################

InfectTimesMedianDFAllBadgerMedian <- median(InfectTimesMedianDFAll[InfectTimesMedianDFAll$host == "Badger",]$infect_time)

InfectTimesMedianDFAllBadgerQuantiles <- quantile(InfectTimesMedianDFAll[InfectTimesMedianDFAll$host == "Badger",]$infect_time, probs=c(0.05, 0.95))

InfectTimesMedianDFAllBovineMedian <- median(InfectTimesMedianDFAll[InfectTimesMedianDFAll$host == "Bovine",]$infect_time)

InfectTimesMedianDFAllBovineQuantiles <- quantile(InfectTimesMedianDFAll[InfectTimesMedianDFAll$host == "Bovine",]$infect_time, probs=c(0.05, 0.95))

##############################
## Relabel bovine to cattle ##
##############################

InfectTimesMedianDFAll <- InfectTimesMedianDFAll %>% 
  mutate(host = str_replace(host, "Bovine", "Cattle"))

######################################################################
## Plot boxplot of infection times coloured by host for all samples ##
######################################################################

infectTimesAllBoxplot <- ggplot() + geom_boxplot(data = InfectTimesMedianDFAll, aes(x=factor(host),y = infect_time, fill=factor(host))) + theme_bw() + ylab("Time infected Before Sampling (Years)") + xlab("Host") + theme(legend.position = "none") + coord_flip()

##############################################################################################
## Plot ridgeplot (Figure 2D) of infection times coloured by host and faceted by cluster id ##
##############################################################################################

infectTimesAllRidgeFacet <- ggplot() + 
  geom_density_ridges(data = transform(InfectTimesMedianDFAll, 
                                       host = factor(host, levels = c("Cattle", "Badger"))), 
                      aes(x = infect_time, 
                          y = as.factor(host),
                          fill = factor(host)),
                      alpha = 0.5,
                      scale = 0.9) +
  scale_fill_manual(values=c("#0000cc", "#ff0000")) +
  theme_bw() + 
  xlab("Time infected Before Sampling (Years)") + 
  ylab("Host") + 
  theme(legend.position = "none") + 
  facet_wrap(~new_network_id, ncol = 3) + 
  theme(strip.text.x = element_text(size = 15, color = "black"),
        strip.background = element_rect(color="black", fill="white", size=1, linetype="solid")) +
  geom_vline(xintercept = InfectTimesMedianDFAllBadgerMedian, colour = "#ff0000", linetype = "dashed", size = 1) + 
  geom_vline(xintercept = InfectTimesMedianDFAllBovineMedian, colour = "#0000cc", linetype = "dashed", size = 1) + 
  theme(axis.title=element_text(size=15)) + 
  theme(axis.text.x = element_text(size = 15)) + 
  theme(axis.text.y = element_text(size = 15))

########################################
## Read in transmission probabilities ##
########################################

################################
## List probability csv files ##
################################

transProbFiles <- list.files(pattern = "*_transmission_probabilities.csv")

###################################
## Extract cluster names as list ##
###################################

transProbFilesNames <- gsub('_transmission_probabilities.csv', '', transProbFiles)

#######################
## Read in csv files ##
#######################

transProbFilesRead <- lapply(transProbFiles, function(x){
  read_csv(x)
})

names(transProbFilesRead) <- transProbFilesNames

##############################################################
## Add pairwise SNP distances to transmission probabilities ##
##############################################################

transProbpwDist <- lapply(transProbFilesRead, function(x){
  left_join(x, M_bovisPWsnpsDecon)
})

#########################################################################
## Collapse transmission probabilities into single df and filter < 0.5 ##
#########################################################################

transProbpwDistFiltered <- bind_rows(transProbpwDist, .id = "network_id")

transProbpwDistFiltered <- transProbpwDistFiltered %>% 
  filter(prob >= 0.5)

##############################################################################
## Add host_link column identifying the host match of each pair of isolates ##
##############################################################################

transProbpwDistFiltered <- transProbpwDistFiltered %>% 
  left_join(allMetadata, by = c("Taxon1" = "lane_id")) %>% 
  left_join(allMetadata, by = c("Taxon2" = "lane_id")) %>% 
  select_if(!grepl("decimal_date", names(.))) %>% 
  mutate(host_link = ifelse(host.x == "Bovine" & host.y == "Bovine", paste("Bovine-Bovine"),
                            ifelse(host.x == "Bovine" & host.y == "Badger", paste("Interspecies"),
                                   ifelse(host.x == "Badger" & host.y == "Bovine", paste("Interspecies"),
                                          paste("Badger-Badger")))))

###############################################################
## Summarize snp distance for likely transmissions (Table 3) ##
###############################################################

transProbpwDistFilteredSummary <- transProbpwDistFiltered %>% 
  group_by(host_link) %>% 
  summarize(median_snp_dist = median(dist, na.rm = TRUE), min_snp_dist = min(dist, na.rm = TRUE), max_snp_dist = max(dist, na.rm = TRUE))

################################
## Read in sampling csv files ##
################################

samplingFiles <- list.files(pattern = "*cases.csv")

###################################
## Extract cluster names as list ##
###################################

samplingFilesNames <- gsub('_incident_cases.csv', '', samplingFiles)

samplingFilesNames <- gsub('network', '', samplingFilesNames)

#######################
## Read in csv files ##
#######################

samplingFilesRead <- lapply(samplingFiles, function(x){
  read_csv(x)
})

names(samplingFilesRead) <- samplingFilesNames

###########################################################################
## Plot sampling distributions for each cluster (Supplementary Figure 4) ##
###########################################################################

samplingFilesReadPlots <- lapply(names(samplingFilesRead), function(x){
  ggplot(samplingFilesRead[[x]]) +
    geom_bar(aes(x = Time, y = Cases, fill = IsSampled), stat="identity") +
    theme_bw() +
    xlab("Date") + 
    ylab("Estimated Number of cases") +
    ggtitle(x) +
    theme(legend.position = "none") +
    theme(axis.title=element_text(size=15)) + 
    theme(axis.text.x = element_text(size = 15)) + 
    theme(axis.text.y = element_text(size = 15)) +
    theme(plot.title = element_text(size=15)) +
    theme(plot.title = element_text(hjust = 0.5))
})

gridExtra::grid.arrange(samplingFilesReadPlots[[1]],
                        samplingFilesReadPlots[[5]],
                        samplingFilesReadPlots[[6]],
                        samplingFilesReadPlots[[7]],
                        samplingFilesReadPlots[[8]],
                        samplingFilesReadPlots[[9]],
                        samplingFilesReadPlots[[10]],
                        samplingFilesReadPlots[[11]],
                        samplingFilesReadPlots[[12]],
                        samplingFilesReadPlots[[2]],
                        samplingFilesReadPlots[[3]],
                        samplingFilesReadPlots[[4]])