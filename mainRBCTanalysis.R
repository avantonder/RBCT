###############################
### Load required libraries ###
###############################

library(tidyverse)
library(rgdal)
library(ggmap)
library(ggrepel)
library(geosphere)
library(patchwork)
library(reshape2)
library(ape)
library(ggtree)
library(lubridate)
library(ade4)
library(ggridges)
library(scatterpie)
library(phytools)

################
### Functions ##
################

#######################################################################
## Function to plot root to tip distances against collection date to ##
## assess temporal signal in dataset                                 ##
#######################################################################

plotRootToTip <- function(datedTree){
  
  # Extract root to tip distances from tree
  
  RootTotipDistances <- diag(vcv.phylo(datedTree))
  
  # Extract tip labels from tree
  
  treeLabels <- datedTree$tip.label
  
  # Split tip dates from tip labels and add to dataframe
  
  treeDF <- colsplit(string = treeLabels, 
                     pattern = "-", 
                     names = c("tip_label", 
                               "date"))
  # Extract tip dates to list
  
  treeDates <- treeDF$date
  
  # Perform linear regression on root to tip distances and tree dates
  
  treeModel <- lm(RootTotipDistances ~ treeDates)
  
  # Calculate correlation between root to tip distances and tree dates
  
  treeCorrelation <- cor.test(treeDates, 
                              RootTotipDistances, 
                              method = "pearson", 
                              conf.level = 0.95)
  
  # Summarise linear regression
  
  modelSummary <- summary(treeModel)
  
  # Calculate x intercept i.e. date of MRCA
  
  xIntercept <- -coef(treeModel)[1]/coef(treeModel)[2]
  
  # Add root to tip distances and tree dates to dataframe for plotting
  
  RootTotipDF <- data.frame(RootTotipDistances,
                            treeDates)
  
  # Create root to tip plot and add variables calculated above
  
  RootTotipPlot <- ggplot(data = RootTotipDF, 
                          aes(treeDates,  
                              RootTotipDistances)) + 
    geom_point(alpha = 0.25, 
               size = 2, 
               colour = "red") + 
    theme_classic() + 
    labs(x = "Year",
         y = "Root to tip distance") +
    geom_smooth(method = 'lm',
                fullrange=T,
                se=T) + 
    ggtitle(paste("Slope: ",
                  formatC(modelSummary$coefficients[2], 
                          format = "e", 
                          digits = 3),
                  "; ", 
                  "TMRCA: ",
                  round(xIntercept, 1),
                  "\n", 
                  "Correlation Coefficient: ",
                  round(treeCorrelation$estimate, 3), 
                  ";", 
                  "R^2: ", 
                  format(modelSummary$r.squared,
                         digits = 3)))
  return(RootTotipPlot)
}

#############################
### Set working directory ###
#############################

ERADbTB_dir <- "~/RBCT"
setwd(ERADbTB_dir)

#######################################################
### Map cattle & badger samples using location data ###
#######################################################

##################################
## Read in APHA metadata file ##
##################################

aphaData <- read_csv("APHA_metadata_090620.csv")

#######################################
## Extract coordinates from above df ##
#######################################

aphaCoords <- aphaData %>%
  select(chmapx,chmapy)

###################################################
## Shortcuts to help convert UK grid to lat-long ##
###################################################

ukgrid <- "+init=epsg:27700"
latlong <- "+init=epsg:4326"

#################################
## Create coordinates variable ##
#################################

coords <- cbind(Easting = as.numeric(as.character(aphaCoords$chmapx)),
                Northing = as.numeric(as.character(aphaCoords$chmapy)))

#######################################
## Create the SpatialPointsDataFrame ##
#######################################

aphaData_SP <- SpatialPointsDataFrame(coords,
                                      data = aphaCoords,
                                      proj4string = CRS("+init=epsg:27700"))

############################################
## Convert from one projection to another ##
############################################

aphaData_SP_LL <- spTransform(aphaData_SP, CRS(latlong))

#######################
## replace Lat, Long ##
#######################

aphaData_SP_LL@data$Long <- coordinates(aphaData_SP_LL)[, 1]
aphaData_SP_LL@data$Lat <- coordinates(aphaData_SP_LL)[, 2]

############################
## Write out as shapefile ##
############################

writeOGR(obj = aphaData_SP_LL, dsn = '.', layer = 'aphaDataPoints', driver = 'ESRI Shapefile')

#######################
## Read in shapefile ##
#######################

aphaShapes <- readOGR(dsn = "aphaDataPoints.shp", stringsAsFactors = F)

#########################################
## Merge aphaShapes df and aphaData df ##
#########################################

aphaData$Long <- aphaShapes@data$Long
aphaData$Lat <- aphaShapes@data$Lat

##########################################################
## Read in new network ids and add them to allShapes df ##
##########################################################

allNewNetworks <- read_csv("ERADbTB_new_networks.csv")

aphaData <- aphaData %>% 
  left_join(allNewNetworks, by = "lane_id")

#########################################
## Define SW UK lat/long for stamenmap ##
#########################################

swUK <- c(left = -6, bottom = 49.5, right = -1, top = 53.5)

##############################################
## Download stamenmap for above coordinates ##
##############################################

stamenmap <- get_stamenmap(swUK, zoom = 8, maptype = "terrain-background", color = "bw", force = TRUE)

####################################################################
## Calculate coordinates for each trial area to place scatterpies ##
####################################################################

TrialAreaCoords <- aphaData %>% 
  group_by(trial_id) %>%
  summarize(median_long = median(Long),
            median_lat = median(Lat)) 

#############################################################
## Read in number of isolates per host for each trial area ##
#############################################################

hostCounts <- read_csv("trial_area_hosts.csv")

TrialAreaCoords <- TrialAreaCoords %>% 
  left_join(hostCounts, by = "trial_id") 

#####################################################################
## Calculate radius for each scatterpie based on number of animals ##
#####################################################################

TrialAreaCoords <- TrialAreaCoords %>%
  mutate(total = Bovines + Badgers) %>% 
  mutate(radius = 0.001 * total)

####################################################################################
## Plot trial area scatterpies showing proportion of isolates by host (Figure 1A) ##
####################################################################################

trialAreaHostScatterPie <- ggmap(stamenmap) + 
  geom_scatterpie(data = TrialAreaCoords,
                  cols = c("Bovines", "Badgers"),
                  aes(x = median_long,
                      y = median_lat,
                      r = radius)) +
  coord_equal() +
  theme(legend.position = "none") +
  geom_scatterpie_legend(radius = TrialAreaCoords$radius, 
                         x = -5.5, 
                         y = 53,
                         n = 3,
                         labeller = function(x) x = TrialAreaCoords$radius[1:5]/0.001) + 
  geom_scatterpie(aes(x=long, y=lat, group=region, r=radius), data=d,
                  cols=LETTERS[1:4], color=NA) + 
  coord_equal()

p + geom_scatterpie_legend(d$radius, x=-140, y=-70)

##############################################################################
## Create facet map of spoligotype clusters sorted by frequency (Figure 1C) ##
##############################################################################

stamenmapSpoligotypeFacet <- ggmap(stamenmap) + 
  geom_point(data = transform(aphaData, 
                              Spoligotype = factor(spoligotype, 
                                                   levels=c("SB0140",
                                                            "SB0263",
                                                            "SB0129",
                                                            "SB0274",
                                                            "SB0957",
                                                            "SB0145"))), 
             aes(x = Long, 
                 y = Lat, 
                 colour = as.factor(Spoligotype), 
                 shape = host), 
             alpha = 0.8) +
  scale_shape_manual(values=c(1, 4)) +
  xlab("Longitude") + 
  ylab("Latitude") + 
  labs(colour = 'Spoligotype',
       shape = 'Host') + 
  facet_wrap(~Spoligotype, 
             ncol=3)

#########################################################
## Create new df to filter out NA new network clusters ##
#########################################################

aphaDataNewNetwork <- aphaData %>% 
  filter(!is.na(new_network_id))

#########################################################################################
## Create facet plot of network clusters sorted by frequency from above df (Figure 2B) ##
#########################################################################################

stamenmapNetworkFacet <- ggmap(stamenmap) + geom_point(data = transform(aphaDataNewNetwork, new_network_id=factor(new_network_id, levels=c(1,2,3,4,5,6,7,8,9,10,12))), 
                                                       aes(x = Long, y = Lat, colour = as.factor(new_network_id), shape = host), alpha = 0.5, fill = NA) + 
  xlab("Longitude") + 
  ylab("Latitude") + 
  labs(colour = 'Transmission cluster', shape = 'Host') + 
  facet_wrap(~new_network_id, ncol=4)

stamenmapNetworkFacet <- ggmap(stamenmap) +
  geom_point(data = aphaDataNewNetwork, 
             aes(x = Long, 
                 y = Lat, 
                 colour = as.factor(new_network_id), 
                 shape = host), 
             alpha = 0.8) + 
  scale_shape_manual(values=c(1, 
                              4)) +
  scale_color_manual(values=c("#ff0000",
                              "#663366",
                              "#ffff00",
                              "#006600",
                              "#9900ff",
                              "#ff6600",
                              "#66ffff",
                              "#aaffaa",
                              "#0000cc",
                              "#ff66cc",
                              "#9999ff",
                              "#996600")) +
  xlab("Longitude") + 
  ylab("Latitude") + 
  labs(colour = 'Transmission cluster', 
       shape = 'Host') + 
  facet_wrap(~new_network_id, 
             ncol=4)

#############################
### Pairwise SNP analyses ###
#############################

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

#################################
## Match sample and trial area ##
#################################

M_bovisPWsnpsDecon$Taxon1_area <- aphaData$trial_id[match(M_bovisPWsnpsDecon$Taxon1, aphaData$lane_id)]
M_bovisPWsnpsDecon$Taxon2_area <- aphaData$trial_id[match(M_bovisPWsnpsDecon$Taxon2, aphaData$lane_id)]

###########################
## Match sample and host ##
###########################

M_bovisPWsnpsDecon$Taxon1_host <- aphaData$host[match(M_bovisPWsnpsDecon$Taxon1, aphaData$lane_id)]
M_bovisPWsnpsDecon$Taxon2_host <- aphaData$host[match(M_bovisPWsnpsDecon$Taxon2, aphaData$lane_id)]

#################################
## Match sample and network_id ##
#################################

M_bovisPWsnpsDecon$Taxon1_network <- aphaData$new_network_id[match(M_bovisPWsnpsDecon$Taxon1, aphaData$lane_id)]
M_bovisPWsnpsDecon$Taxon2_network <- aphaData$new_network_id[match(M_bovisPWsnpsDecon$Taxon2, aphaData$lane_id)]

##############################################################
## Create new df with columns for area_match and host_match ##
## to help with separating and filtering data for plots     ##
##############################################################

M_bovisPWsnpsDeconFiltered <- M_bovisPWsnpsDecon %>%
  filter(Taxon1_area != "Unknown" & Taxon2_area != "Unknown") %>% 
  mutate(area_match = ifelse(Taxon1_area == Taxon2_area, paste("same_trial_area"), paste("different_trial_area"))) %>% 
  mutate(host_match = ifelse(Taxon1_host == Taxon2_host, paste("same_host"), paste("different_host"))) %>% 
  mutate(network_match = ifelse(Taxon1_network == Taxon2_network, paste("Same Transmission Cluster"), paste("Different Transmission Cluster"))) %>%
  mutate(host_link = ifelse(Taxon1_host == "Bovine" & Taxon2_host == "Bovine", paste("Cattle-Cattle"),
                            ifelse(Taxon1_host == "Bovine" & Taxon2_host == "Badger", paste("Cattle-Badger"),
                                   ifelse(Taxon1_host == "Badger" & Taxon2_host == "Bovine", paste("Cattle-Badger"),
                                          paste("Badger-Badger")))))

#################################################################
## Create new df filtering out samples that aren't in the same ##
## trial area                                                  ##
#################################################################

M_bovisPWsnpsDeconFilteredClusters <- M_bovisPWsnpsDeconFiltered %>% 
  filter(!is.na(network_match)) 

#######################################################################
## Plot all pairwise distances as histogram (Supplementary Figure 1) ##
#######################################################################

M_bovisHist1 <- ggplot(data = M_bovisPWsnpsDeconFilteredClusters, 
                       aes(x = dist, fill = network_match)) + 
  geom_histogram(bins = 200, 
                 alpha = 0.8, 
                 position = "identity") + 
  theme_bw()  + 
  xlab("Pairwise SNP distance") + 
  ylab("Frequency") + 
  facet_wrap(~host_link, ncol=1) + 
  scale_fill_discrete(name = "Transmission cluster") + 
  theme(axis.title=element_text(size=15)) + 
  theme(axis.text.x = element_text(size = 15)) + 
  theme(axis.text.y = element_text(size = 15)) +
  theme(strip.text.x = element_text(size = 15, color = "black"),
        strip.background = element_rect(color="black", fill="white", size=1, linetype="solid"))

##########################################
### Geographical localization analyses ###
##########################################

################################################################
## Create a copy of deconvoluted pairwise SNP distance matrix ##
################################################################

ERADbTBpwSNPsGeoDists <- M_bovisPWsnpsDecon

###########################################################################
## Create a new dataframe with lane_id, Spoligotype, Longitude, Latitude ##
###########################################################################

allShapesGeoDists <- aphaData %>% 
  select(lane_id, spoligotype, Long, Lat,new_network_id, trial_id)

#####################################
## Add spoligotypes for each taxon ##
#####################################

ERADbTBpwSNPsGeoDists <- ERADbTBpwSNPsGeoDists %>% 
  left_join(allShapesGeoDists, by = c("Taxon1" = "lane_id")) %>% 
  left_join(allShapesGeoDists, by = c("Taxon2" = "lane_id"))

######################################################
## Calculate pw geographical distance for each pair ##
######################################################

ERADbTBpwSNPsGeoDists <- ERADbTBpwSNPsGeoDists %>% 
  mutate(geo_dist_km = distHaversine(cbind(Long.x, Lat.x), cbind(Long.y, Lat.y))/1000) %>% 
  mutate(geo_dist_miles = 0.621371 * geo_dist_km)

#################################################################
## Create new dataframe for geographical localization analyses ##
#################################################################

ERADbTBpwSNPsGeoDistsLoc <- ERADbTBpwSNPsGeoDists %>% 
  select(Taxon1, Taxon2, spoligotype.x, spoligotype.y, dist, geo_dist_km) %>% 
  mutate(spoligo_match = ifelse(spoligotype.x == spoligotype.y, paste("same_spoligo"), paste("different_spoligo"))) %>% 
  filter(spoligo_match == "same_spoligo") %>% 
  group_by(spoligotype.x) %>%
  filter(n() >= 490) %>%
  ungroup()

################################################################################
## Plot Pairwise SNP distribution for most prevalent spoligotypes (Figure 1D) ##
################################################################################

ERADbTBpwSNPsGeoDistsLocDistRidges <- ggplot(data = transform(ERADbTBpwSNPsGeoDistsLoc, 
                                                              Spoligotype = factor(spoligotype.x,
                                                                                   levels=c("SB0145",
                                                                                            "SB0957",
                                                                                            "SB0274",
                                                                                            "SB0129",
                                                                                            "SB0263",
                                                                                            "SB0140"))), 
                                             aes(x = dist, 
                                                 y = Spoligotype, 
                                                 fill = Spoligotype)) +
  geom_density_ridges(scale = 0.9) +
  scale_fill_manual(values=rev(c("#0000cc", 
                                 "#ff0000", 
                                 "#ff6600", 
                                 "#006600", 
                                 "#994d00", 
                                 "#ff66cc"))) +
  theme_bw() +
  labs(x = "Pairwise SNP distance") +
  theme(axis.title=element_text(size=15)) + 
  theme(axis.text.x = element_text(size = 15)) + 
  theme(axis.text.y = element_text(size = 15)) +
  theme(legend.position = "none")

##########################################################################################################
## Plot comparison of pairwise SNP and geographic distances for most prevalent spoligotypes (Figure 1E) ##
##########################################################################################################

ERADbTBpwSNPsGeoDistsLocGeoDistDistPlot <- ggplot(data = transform(ERADbTBpwSNPsGeoDistsLoc, Spoligotype=factor(spoligotype.x, levels=c("SB0140","SB0263","SB0129","SB0274","SB0957","SB0145"))), aes(x = dist, y = geo_dist_km)) +
  geom_point(aes(colour = Spoligotype), size = 0.4) +
  scale_color_manual(values=c("#0000cc", "#ff0000", "#ff6600", "#006600", "#994d00", "#ff66cc")) +
  theme_bw() +
  xlab("Pairwise SNP distance") + 
  ylab("Geographic Distance (km)") + 
  facet_wrap(~Spoligotype, ncol=3) + 
  theme(axis.title=element_text(size=15)) + 
  theme(axis.text.x = element_text(size = 15)) + 
  theme(axis.text.y = element_text(size = 15)) +
  theme(legend.position = "none") + 
  theme(strip.text.x = element_text(size = 15, color = "black"),
        strip.background = element_rect(color="black", fill="white", size=1, linetype="solid"))

#######################################
### Network definition and analyses ###
#######################################

###############################################################
### Create network of all samples using SNP threshold of 15 ###
###############################################################

################################################
## Create new df containing nodes for network ##
################################################

allShapes_nodes_15_snps <- as_tibble(allShapes) %>%
  mutate(id = row_number()) %>% 
  select(id, lane_id, chtype, trial_area, fastbaps, Long, Lat)

#############################################################
## Extract pairwise SNP comparisons < 16 (edges of network) ##
#############################################################

M_bovisPWsnpsDeconFiltered_edges_15_snp <- M_bovisPWsnpsDeconFiltered %>%
  filter(dist < 16) %>% 
  select(Taxon1, Taxon2, dist)

M_bovisPWsnpsDeconFiltered_edges_15_snp <- M_bovisPWsnpsDeconFiltered_edges_15_snp %>%
  left_join(allShapes_nodes_15_snps[,c(1,2)], by = c("Taxon1" = "lane_id")) %>% 
  dplyr::rename(from.id = id)

M_bovisPWsnpsDeconFiltered_edges_15_snp <- M_bovisPWsnpsDeconFiltered_edges_15_snp %>% 
  left_join(allShapes_nodes_15_snps[,c(1,2)], by = c("Taxon2" = "lane_id")) %>% 
  select(from.id, to.id = id, dist)

###########################
## Create network object ##
###########################

all_15_snp_routes <- tbl_graph(nodes = allShapes_nodes_15_snps, edges = M_bovisPWsnpsDeconFiltered_edges_15_snp, directed = TRUE)

all_15_snp_routes_components <- components(all_15_snp_routes)

allShapes_nodes_15_snps_networked <- allShapes_nodes_15_snps %>% 
  mutate(network_id = all_15_snp_routes_components$membership) 

M_bovisPWsnpsDeconFiltered_networked <- M_bovisPWsnpsDeconFiltered

M_bovisPWsnpsDeconFiltered_networked$Taxon1_network <- allShapes_nodes_15_snps_networked$network_id[match(M_bovisPWsnpsDeconFiltered_networked$Taxon1, allShapes_nodes_15_snps_networked$lane_id)]
M_bovisPWsnpsDeconFiltered_networked$Taxon2_network <- allShapes_nodes_15_snps_networked$network_id[match(M_bovisPWsnpsDeconFiltered_networked$Taxon2, allShapes_nodes_15_snps_networked$lane_id)]

#############################################
### Root to tip analyses for all networks ###
#############################################

########################################################
## Create new dataframe with dates and new tip labels ##
########################################################

sampleDates <- aphaData %>% 
  select(lane_id, date_taken) %>% 
  mutate(decimalDate = decimal_date(dmy(date_taken))) %>% 
  mutate(newTipLabel = paste(lane_id, 
                             decimalDate, 
                             sep = '-'))

############################
## List iqtree tree files ##
############################

allNetworkTreeFiles <- list.files(pattern = "*_masked_snps.aln.treefile")

###################################
## Extract network names as list ##
###################################

allNetworkTreeFilesNames <- gsub('_130320_masked_snps.aln.treefile', '', allNetworkTreeFiles)

allNetworkTreeFilesNames <- gsub('network', '', allNetworkTreeFilesNames)

###################################################
## Read in tree files. Add network names to list ##
###################################################

allNetworkTreeFilesRead <- lapply(allNetworkTreeFiles, function(x){
  read.newick(x)
})

names(allNetworkTreeFilesRead) <- allNetworkTreeFilesNames

#############################################################################
## Add collection dates to tip labels of trees and midpoint root the trees ##
#############################################################################

allNetworkTreeFilesDated <- allNetworkTreeFilesRead

for (i in 1:length(allNetworkTreeFilesRead)){
  
  allNetworkTreeFilesDated[[i]]$tip.label <- sampleDates$newTipLabel[match(allNetworkTreeFilesDated[[i]]$tip.label,
                                                                           sampleDates$lane_id)]
}


allNetworkTreeFilesRooted <- lapply(allNetworkTreeFilesDated, function(x){
  midpoint.root(x)
})

############################################################################
## Plot root to distance against collection date (Supplementary Figure 2) ##
############################################################################

allNetworkR2TPlot <- lapply(allNetworkTreeFilesRooted, function(x){
  plotRootToTip(x)
})

gridExtra::grid.arrange(grobs = allNetworkR2TPlot)

gridExtra::grid.arrange(allNetworkR2TPlot[[1]],
                        allNetworkR2TPlot[[5]],
                        allNetworkR2TPlot[[6]],
                        allNetworkR2TPlot[[7]],
                        allNetworkR2TPlot[[8]],
                        allNetworkR2TPlot[[9]],
                        allNetworkR2TPlot[[10]],
                        allNetworkR2TPlot[[11]],
                        allNetworkR2TPlot[[12]],
                        allNetworkR2TPlot[[2]],
                        allNetworkR2TPlot[[3]],
                        allNetworkR2TPlot[[4]])

##############################################################################
### Plot sampling times and BEAST MRCA dates for each transmission cluster ###
##############################################################################

######################################
## Read in BEAST MRCA dates and CIs ##
######################################

networkBeastDates <- read_csv("BEAST_dates.csv")

networkBeastDates <- networkBeastDates %>% 
  mutate(new_mrca = as.Date(date_decimal(mrca))) %>% 
  mutate(new_lower_CI = as.Date(date_decimal(lower_CI))) %>%
  mutate(new_higher_CI = as.Date(date_decimal(higher_CI))) %>% 
  mutate(date_range = paste(round(mrca, digits = 0),
                            " [",
                            round(lower_CI, digits = 0),
                            "-",
                            round(higher_CI, digits = 0),
                            "]",
                            sep = ""))

################################################################################
## Create new df from aphaDataNewNetwork. Convert collection dates to decimal ##
################################################################################

aphaDataNewNetworkDates <- aphaDataNewNetwork %>% 
  select(lane_id, host, trial_id,date_taken, new_network_id)

#######################################################
## Plot network MRCAs and sampling dates (Figure 2C) ##
#######################################################

aphaDataNewNetworkDatesRidges <- ggplot() + 
  geom_density_ridges(data = subset(transform(aphaDataNewNetworkDates, 
                                              new_network_id = factor(new_network_id, 
                                                                      levels=c(12,11,10,9,8,7,6,5,4,3,2,1))), host == "Badger"), 
                      aes(x = dmy(date_taken), 
                          y = as.factor(new_network_id)),
                      alpha = 0.5,
                      fill = "#ff0000",
                      scale = 0.9) +
  geom_density_ridges(data = subset(transform(aphaDataNewNetworkDates, 
                                              new_network_id = factor(new_network_id, 
                                                                      levels=c(12,11,10,9,8,7,6,5,4,3,2,1))), host == "Bovine"), 
                      aes(x = dmy(date_taken), 
                          y = as.factor(new_network_id)),
                      alpha = 0.5,
                      fill = "#0000cc",
                      scale = 0.9) +
  geom_point(data = networkBeastDates, 
             aes(x = new_mrca, 
                 y = as.factor(network)),
             size = 4) +
  geom_segment(data = networkBeastDates, 
               aes(x = new_lower_CI, 
                   y = as.factor(network),
                   xend = new_higher_CI,
                   yend = as.factor(network)),
               size = 2) +
  geom_text(data = networkBeastDates, 
            aes(x = new_mrca, 
                y = as.factor(network),
                label = date_range),
            position = position_nudge(y = 0.2), 
            size = 4) +
  theme_bw() + 
  xlab("Sampling date") + 
  ylab("Transmission cluster") + 
  geom_vline(xintercept = ymd("2001-02-01"), 
             colour = "black", 
             linetype = "dashed") +
  geom_vline(xintercept = ymd("2001-12-31"), 
             colour = "black", 
             linetype = "dashed") + 
  theme(axis.title=element_text(size=15)) + 
  theme(axis.text.x = element_text(size = 15)) + 
  theme(axis.text.y = element_text(size = 15)) 