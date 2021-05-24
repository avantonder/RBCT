##################################################################################
## Script for parsing and summarising data extracted from CTS and Sam databases ##
##################################################################################

################################
## Read in required libraries ##
################################

library(tidyverse)
library(lubridate)
library(ggtree)
library(treeio)
library(ggnewscale)
library(ggsn)
library(RColorBrewer)
library(patchwork)
library(igraph)
library(tidygraph)
library(ggraph)
library(mapdata)
library(maps)
library(ggrepel)

#############################
## Read in cattle contacts ##
#############################

cattle_contacts <- read_csv("Raw_data/Generated_data/cattle_contacts.csv")

###############################
## Read in animal identities ##
###############################

animal_identities <- read_csv("Raw_data/Generated_data/animal_identities.csv")

APHA_metadata <- animal_identities %>% 
  dplyr::rename(lane_id = LaneID) %>% 
  select(lane_id, everything())

#################################################
## Combine animal identies and final locations ##
#################################################

APHA_metadata_finallocations <- animal_identities %>% 
  dplyr::rename(lane_id = LaneID) %>% 
  select(lane_id, everything()) %>% 
  left_join(final_cattle_locations[,c("onlocationkey", "lane_id")], by = "lane_id") %>% 
  dplyr::rename(locationkey = onlocationkey)

########################
## Filter out badgers ##
########################

bovine_identities <- animal_identities %>% 
  filter(host == "Bovine")

###########################################
## Import movements of sequenced animals ##
###########################################

seq_animal_moves <- read_csv("Raw_data/Generated_data/seq_animal_moves.csv", 
                             col_types = c(importcountry = "c", importdate = "T"), na = c("", "NA"))

######################################################
## Import 7 day movement data for sequenced animals ##
######################################################

seq_animal_moves_continuous_7_days <- read_csv("Raw_data/Generated_data/seq_animal_moves_continuous_7_days.csv",
                                               col_types = c(importcountry = "c", importdate = "T"))

seq_animal_moves_continuous_7_days <- seq_animal_moves_continuous_7_days %>% 
  mutate(animalid = as.character(animalid))

#####################################
## Import movements of all animals ##
#####################################

all_animal_moves_continuous_7_days <- read_csv("Raw_data/Generated_data/All_animal/all_animal_moves_continuous_7_days.csv",
                                               col_types = c(standard_animal_id = "c"))

###########################################
## Import final locations of all animals ##
###########################################

animal_final_locations <- read_csv("Raw_data/Generated_data/animal_final_locations.csv")

animal_final_locations_plot <- animal_final_locations %>% 
  dplyr::rename(trial_label = trial_id)

######################################
## Import final locations of cattle ##
######################################

final_cattle_locations <- read_csv("Raw_data/Generated_data/final_cattle_locations.csv")

################################################################################
## Import fuzzed geographical coordinates for animals not located in database ##
################################################################################

fuzzed_missing_LatLong <- read_csv("Raw_data/Generated_data/fuzzed_missing_LatLong.csv")

missing_animal_plot <- fuzzed_missing_LatLong %>% 
  dplyr::rename(trial_label = trial_id) 

############################################################
## Import fuzzed geographical coordinates for all animals ##
############################################################

fuzzed_Lat_Long_allanimal <- read_csv("Raw_data/Generated_data/All_animal/fuzzed_Lat_Long_allanimal.csv")

#######################################################
## Import geographical coordinates for all locations ##
#######################################################

location_Lat_Longs <- read_csv("Raw_data/Generated_data/location_Lat_Longs.csv")

###########################################################
## Import metadata from BIGSdb for all sequenced animals ##
###########################################################

all_sequenced_animals <- read_csv("Raw_data/ERADbTB_metadata_200720_UPDATED_METADATA.csv") %>% 
  select(`Supplier name`, animal_id, LaneID, host, trial_id, date_taken, Year, chmapx, chmapy, cphh, chcphh)

###################################################
## Extract badger locations from BIGSdb metadata ##
###################################################

badger_locations <- all_sequenced_animals %>% 
  select(-cphh, -chcphh) %>% 
  filter(host == "Badger") %>% 
  rename(eartag = `Supplier name`) %>% 
  left_join(complete_new_networks, by = c("LaneID" = "lane_id")) %>% 
  dplyr::rename(animalid = animal_id) %>% 
  mutate(animalid = if_else(is.na(animalid), str_remove_all(eartag, pattern = "[-/ ]"), animalid))

####################################################
## Extract badger identities from BIGSdb metadata ##
####################################################

badger_identities <- badger_locations %>% 
  mutate(corrected_standardeartag = NA) %>% 
  select(eartag, corrected_standardeartag, animalid, LaneID, host, trial_id, new_network_id, date_taken)

#############################################
## Import transmission cluster information ##
#############################################

complete_new_networks <- read_csv("Raw_data/Generated_data/complete_new_networks.csv")

#########################################################
## Import inferred dates of infection for each cluster ##
#########################################################

infectiontimesfiles <- list.files(path = "Raw_data/Inferred_infection_times/")

network_infectiontimes = data.frame(sample = NA, mean = NA, sdev = NA, min = NA, max = NA, median = NA, UQ = NA, LQ = NA, network_id = NA, animalid = NA)

for (file in infectiontimesfiles) {
  
  file_with_path <- paste("Raw_data/Inferred_infection_times/", file, sep = "")
  
  file_data <- read_csv(file_with_path)
  
  assign(paste(str_remove(file, pattern = "_beast_infection_times_1.csv"), "_infectiontimes", sep = ""), file_data)
  
  summary <- data.frame(sample = colnames(file_data), mean = NA, sdev = NA, min = NA, max = NA, median = NA, UQ = NA, LQ = NA, network_id = NA) # removed ', mode = NA'
  
  summary$mean <- sapply(file_data, mean) 
  summary$sdev <- sapply(file_data, sd)
  summary$min <- sapply(file_data, min)
  summary$max <- sapply(file_data, max)
  summary$median <- sapply(file_data, median)
  summary$UQ <- sapply(file_data, quantile, probs = 0.75)
  summary$LQ <- sapply(file_data, quantile, probs = 0.25)
  #summary$mode <- sapply(network11_infectiontimes, Modes) - problem of distributions could be bimodel etc
  
  summary_V2 <- summary %>% 
    left_join(animal_identities[,c("animalid", "LaneID")], by = c("sample" = "LaneID")) %>% 
    mutate(network_id = capture.output(cat(str_remove(str_remove(file, pattern = "_beast_infection_times_1.csv"), pattern = "network"), sep = "")))
  
  assign(paste(str_remove(file, pattern = "_beast_infection_times_1.csv"), "_infection_summary", sep = ""), summary_V2)
  
  network_infectiontimes <- rbind(network_infectiontimes, summary_V2) %>% 
    filter(!is.na(mean))
  
  rm(file_with_path, file_data, summary, summary_V2)
  
}

#######################################
## Import SNP trees for each cluster ##
#######################################

for (file in list.files(path = "Raw_data/snp_trees/")){
  
  file_with_path <- paste("Raw_data/snp_trees/", file, sep = "")
  
  tree <- ggtree::read.tree(file_with_path)
  
  file_name <- paste(str_remove_all(substr(file, start = 41, stop = 49), pattern = "[_]"), "_tree", sep = "")
  
  assign(file_name, tree)
  
  
  rm(file_with_path, tree, file_name)
  
}

#################################################
## Read in transmission nodes for each cluster ##
#################################################

transNodesFiles <- list.files(path = "Raw_data/Transmission_networks/", pattern = "*Nodes.csv")

###################################
## Extract network names as list ##
###################################

transNodesFilesNames <- gsub('TransmissionNodes.csv', '', transNodesFiles)

#######################
## Read in csv files ##
#######################

transNodesFilesRead <- lapply(transNodesFiles, function(x){
  read_csv(paste("Raw_data/Transmission_networks/", x, sep = ""))
})

names(transNodesFilesRead) <- transNodesFilesNames

#################################################
## Read in transmission edges for each cluster ##
#################################################

transEdgesFiles <- list.files(path = "Raw_data/Transmission_networks/", pattern = "*Edges.csv")

###################################
## Extract network names as list ##
###################################

transEdgesFilesNames <- gsub('TransmissionEdges.csv', '', transEdgesFiles)

#######################
## Read in csv files ##
#######################

transEdgesFilesRead <- lapply(transEdgesFiles, function(x){
  read_csv(paste("Raw_data/Transmission_networks/", x, sep = ""))
})

names(transEdgesFilesRead) <- transEdgesFilesNames

################################################
## Add metadata to transmission cluster nodes ##
################################################

transNodes <- NULL

for (i in 1:length(transNodesFilesRead)){
  transNodes[[i]] <- transNodesFilesRead[[i]] %>% 
    dplyr::rename(metadata_trial_id = trial_id) %>% 
    left_join(animal_final_locations[,c("lane_id", "animalid", "onlocationkey", "trial_id")], by = c("label" = "lane_id")) %>% 
    left_join(fuzzed_missing_LatLong[,c("lane_id", "locationkey")], by = c("label" = "lane_id")) %>% 
    mutate(onlocationkey = ifelse(is.na(onlocationkey), locationkey, onlocationkey)) %>% 
    select(-locationkey) %>% 
    mutate(trial_id = ifelse(is.na(trial_id), metadata_trial_id, trial_id))
}

names(transNodes) <- transEdgesFilesNames

##################################
## Create transmission networks ##
##################################

transNetworksNames <- transEdgesFilesNames

TransmissionNetworks <- NULL
for (i in 1:length(transNodes)){
  TransmissionNetworks[[i]] <- graph_from_data_frame(transEdgesFilesRead[[i]], vertices = transNodes[[i]]) %>% as_tbl_graph()
}

names(TransmissionNetworks) <- transNetworksNames

################################################
## Assign colours to trial areas for plotting ##
################################################

trial_id_colours <- c("A3" = "#e41a1c", "B2" = "#1f78b4", "C3" = "#1b9e77", "D3" = "#d95f02", "E3" = "#7570b3", "F1" = "#e7298a", "G2" = "#66a61e",
                      "H2" = "#e6ab02", "I2" = "#a6761d", "J1" = "#666666", "other" = "black")

##########################
## Create UK map object ##
##########################

UK <- map_data(map = "worldHires",
               region = c("UK:Great Britain", "Isle of Wight", "Wales:Anglesey", "Isle of Man", "UK:Isle of Sheppey", "UK:Saint Mary's", "UK:Scotland:Island of Arran",
                          "UK:Scotland:Colonsay", "UK:Scotland:Island of Mull", "UK:Scotland:North Uist", "UK:Scotland:South Uist", "UK:Scotland:Islay",
                          "UK:Scotland:Jura", "UK:Scotland:Barra", "UK:Scotland:Benbecula", "UK:Scotland:Tiree", "UK:Scotland:Coll", "UK:Scotland:Island of Skye", "UK:Scotland:Ruhm", "UK:Scotland:Isle of Lewis"))

######################################
## Assign shapes to different hosts ##
######################################

shapes_animals <- c(4, 16)

names(shapes_animals) <- c("Badger", "Bovine")

######################################
## Assign colours to life histories ##
######################################

death_IQR_colours <- c("Death" = "red", "IQR of infection date" = "green", "Birth" = "blue")

################
## Recurrence ##
################

#################################
#### Cluster 6 - location 6A ####
#################################

#############################################
## Extract cattle contacts for location 6A ##
#############################################

cattle_contacts_6A <- cattle_contacts %>% 
  filter(onlocationkey == "6A") %>% 
  mutate(epi_group = NA) %>% 
  mutate(epi_group = ifelse(deathdate < 2001.12, 1, epi_group)) %>% 
  mutate(epi_group = ifelse(deathdate > 2001.12 & deathdate < 2001.4, 2, epi_group)) %>% 
  mutate(epi_group = ifelse(deathdate > 2001.4 & deathdate < 2003.6, 3, epi_group)) %>%
  mutate(epi_group = ifelse(deathdate > 2003.6 & deathdate < 2003.9, 4, epi_group)) %>% 
  mutate(epi_group = ifelse(deathdate > 2003.9 & deathdate < 2004.7, 5, epi_group)) %>% 
  mutate(epi_group = ifelse(deathdate > 2004.7, 6, epi_group)) %>% 
  arrange(epi_group) %>% 
  mutate(rownumber = row_number()) %>% 
  mutate(animalid = as.character(animalid)) %>% 
  mutate(type = "Death")

################################################################
## Extract infection times from cluster 6 infection dataframe ##
################################################################

network6_infection_summary_6A <- network6_infection_summary %>% 
  filter(animalid %in% cattle_contacts_6A$animalid) %>% 
  left_join(cattle_contacts_6A[,c("animalid", "rownumber")], by = "animalid") %>% 
  mutate(type = "IQR of infection date")

###################################################
## Extract birth dates for cattle in location 6A ##
###################################################

cattle_birthdate_6A <- cattle_contacts %>% 
  filter(animalid %in% cattle_contacts_6A$animalid) %>% 
  filter(birth == TRUE) %>% 
  mutate(animalid = as.character(animalid)) %>% 
  left_join(cattle_contacts_6A[,c("animalid", "rownumber")], by = "animalid") %>% 
  mutate(type = "Birth")

###########################################
## Define herd breakdowns at location 6A ##
###########################################

breakdowns2 <- data.frame(start_Date = c(2001.022, 2003.499, 2005.734), end_Date = c(2001.901, 2005.515, 2006.677))

#####################################
## Assign slaughter group colours  ##
#####################################

slaughter_group_colours <- c("1" = "#d6604d", "2" = "#f4a582", "3" = "#d1e5f0", "4" = "#92c5de", "5" = "#4393c3", "6" = "#2166ac")

#############################################################################
## Plot life histories for animals involved in herd breakdowns (Figure 3A) ##
#############################################################################

plot1_6A <- ggplot() +
  geom_segment(data = network6B_infection_summary_6A, aes(x = LQ, xend = LQ, y = as.factor(rownumber), yend = as.factor(rownumber), colour = NULL)) +
  geom_segment(data = network6B_infection_summary_6A, aes(x = LQ, xend = LQ, y = rownumber - 0.3, yend = rownumber + 0.3, colour = type), size = 1)+
  geom_segment(data = network6B_infection_summary_6A, aes(x = UQ, xend = UQ, y = rownumber - 0.3, yend = rownumber + 0.3, colour = type), size = 1)+
  geom_segment(data = cattle_contacts_6A, aes(x = movementdate, xend = leaving_date, y = rownumber, yend = rownumber))+
  geom_segment(data = cattle_contacts_6A, aes(x = deathdate, xend = deathdate, y = rownumber-0.3, yend = rownumber + 0.3, colour = type), size = 1) +
  geom_segment(data = cattle_birthdate_6A, aes(x= movementdate, xend = movementdate, y = rownumber-0.3, yend = rownumber+0.3, colour = type), size = 1)+
  geom_rect(data = breakdowns2, aes(xmin = start_Date, xmax = end_Date, ymin = 0, ymax = 27), fill = "#666666", alpha = 0.2) +
  geom_tile(data = cattle_contacts_6A, aes(x = 2007.2, y = rownumber, fill = as.factor(epi_group), width = 0.8), show.legend = FALSE)+
  scale_fill_manual(values = slaughter_group_colours)+
  scale_colour_manual(values = death_IQR_colours, labels = c("Birth", "Death", "IQR of\ninfection date")) +
  theme_bw() +
  xlab("Date")+
  ylab("Animal ID")+
  labs(fill = "Slaughter group", colour = "Event") +
  scale_x_continuous(breaks = c(1998,2000,2002,2004,2006), labels = c(1998,2000,2002,2004,2006))+
  guides(colour = guide_legend(order = 1)) + 
  theme(axis.title=element_text(size=15)) + 
  theme(axis.text.x = element_text(size = 15)) +
  theme(axis.text.y = element_text(size = 15)) +
  theme(legend.position = "none")

#################################################################
## Prepare dataframe for adding metadata to cluster 6 SNP tree ##
#################################################################

tree_lab_6A <- cattle_contacts_6A %>% 
  left_join(bovine_identities[,c("LaneID", "animalid")], by = "animalid") %>% 
  select(LaneID, animalid, epi_group) %>% 
  dplyr::rename(lane_id = LaneID) %>% 
  dplyr::rename(epi_label = epi_group)

tree_lab_6ADF <- as.data.frame(tree_lab_6A[,3])

rownames(tree_lab_6ADF) = tree_lab_6A$lane_id

tree_lab_6ADF$epi_label <- as.factor(tree_lab_6ADF$epi_label)

APHA_metadata_psudonym_6A <- animal_identities %>% 
  left_join(cattle_contacts_6A[,c("animalid", "rownumber")], by = "animalid") %>% 
  dplyr::rename(lane_id = LaneID) %>% 
  select(lane_id, everything())

#########################################
## Plot cluster 6 SNP tree (Figure 3B) ##
#########################################

network6_plotted <- ggtree(network6_tree)%<+% APHA_metadata_psudonym_6A + 
  theme(legend.position="right") + 
  geom_tippoint(aes(color = host), size=1.5) +
  scale_color_manual(values = c("#ff0000", "#0000cc")) +
  labs(color = "Host animal")+
  geom_tiplab(aes(subset = animalid %in% tree_lab_6A$animalid, label = rownumber), align = T, size = 3, hjust = -0.2)+
  theme_tree2() +
  xlab("SNP distance") +
  theme(text = element_text(size = 15), axis.text.x = element_text(colour = "black")) + 
  xlim_tree(16)

###########################################
## Add slaughter group as metadata strip ##
###########################################

network6_plottedMeta <- gheatmap(network6_plotted, tree_lab_6ADF, width = 0.2, colnames = F, color = NULL) +
  scale_fill_manual(values = slaughter_group_colours, na.translate=FALSE, name = "Slaughter group") +
  labs(color = "Host animal", fill = "Slaughter group") +
  theme(legend.position = "none")

######################
### Superspreading ###
######################

################
## Cluster 12 ##
################

#####################################################
## Extract movement data for animals in cluster 12 ##
#####################################################

cattle_movements_12 <- cattle_contacts %>% 
  dplyr::filter(new_network_id == "12") %>% 
  mutate(animalid = as.character(animalid)) %>% 
  left_join(bovine_identities[,c("LaneID", "animalid")], by = "animalid") %>% 
  select(LaneID, everything()) %>% 
  dplyr::rename(lane_id = LaneID) %>%
  select(lane_id, animalid, onlocationkey, movementdate, leaving_date) %>% 
  arrange(animalid, leaving_date) %>% # from here downwards am trying to create a way of shifting the y position as animals move to different locations
  group_by(animalid) %>% 
  mutate(y_pos = row_number()/10 - 0.1) %>% 
  ungroup()

###############################################
## Identify largest movement for each animal ##
###############################################

max_shift_12 <- cattle_movements_12 %>% 
  select(lane_id, y_pos) %>% 
  arrange(lane_id, desc(y_pos)) %>% 
  group_by(lane_id) %>% 
  mutate(selector = row_number()) %>% 
  ungroup() %>% 
  filter(selector == 1) %>% 
  select(-selector)

####################################################
## Add in coordinates for animals not in database ##
####################################################

missing_animal_12 <- fuzzed_missing_LatLong %>% 
  filter(new_network_id == "12")%>% 
  mutate(trial_label = trial_id) %>% 
  dplyr::rename(onlocationkey = locationkey)

#############################################
## Extract death dates for missing animals ##
#############################################

missing_animal12_deathdate <- missing_animal_12 %>% 
  select(lane_id, deathdate) %>% 
  mutate(animalid = NA, colour_label = "Death", y_pos = 0) %>% 
  select(lane_id, animalid, deathdate, y_pos, colour_label)

#########################################
## Extract death dates for all animals ##
#########################################

death_date_12 <- cattle_contacts %>% 
  dplyr::filter(new_network_id == "12") %>%
  mutate(animalid = as.character(animalid)) %>%
  select(animalid, deathdate) %>% 
  unique() %>% 
  left_join(bovine_identities[,c("LaneID", "animalid")], by = "animalid") %>% 
  select(LaneID, everything()) %>% 
  dplyr::rename(lane_id = LaneID) %>% 
  left_join(max_shift_12, by = "lane_id") %>% 
  mutate(colour_label = "Death") %>% 
  rbind(missing_animal12_deathdate)

#####################################################
## Extract inferred infection times for cluster 12 ##
#####################################################

inferred_infection_time_12 <- network12_infection_summary[,c("sample", "animalid", "median", "LQ", "UQ", "network_id")] %>% 
  filter(network_id == "12") %>% 
  dplyr::rename(lane_id = sample) %>% 
  left_join(max_shift_12, by = "lane_id")%>% 
  mutate(y_pos = ifelse(is.na(y_pos), 0, y_pos)) %>% 
  mutate(colour_label = "IQR of infection date")

############################################
## Extract badger metadata for cluster 12 ##
############################################

badgers_12 <- badger_identities %>% 
  filter(new_network_id == "12") %>% 
  dplyr::rename(lane_id = LaneID) %>% 
  dplyr::rename(badger_date = date_taken) %>% 
  select(lane_id, host, badger_date, trial_id) %>% 
  mutate(colour_label = "Death") %>% 
  mutate(badger_date = decimal_date(as.Date(badger_date, format = "%Y-%m-%d")))

################################################
## Specify colour palette for location colours #
################################################

myColors <- brewer.pal(n = 44, "YlGnBu")

myColorsnames <- cattle_movements_12 %>% 
  select(lane_id, onlocationkey) %>%
  rbind(missing_animal_12[,c("lane_id", "onlocationkey")]) %>% 
  unique() %>% 
  select(onlocationkey) %>% 
  unique()

myColorsnames <- myColorsnames[["onlocationkey"]]

names(myColors) <- myColorsnames

###########################################
## Assign network hub names for plotting ##
###########################################

hubs <- c("1", "2")

#######################################################
## Extract final locations for animals in cluster 12 ##
#######################################################

network12Samples <- animal_final_locations_plot %>% 
  filter(new_network_id == "12")

#################################################
## Extract final locations for missing animals ##
#################################################

network12Missing <- missing_animal_plot %>% 
  filter(new_network_id == "12")

##################################################################
## Prepare dataframe for adding metadata to cluster 12 SNP tree ##
##################################################################

network12SamplesDF <- as.data.frame(network12Samples[,6])

network12MissingDF <- as.data.frame(network12Missing[,11])

rownames(network12SamplesDF) = network12Samples$lane_id

rownames(network12MissingDF) = network12Missing$lane_id

network12All <- rbind(network12SamplesDF, network12MissingDF)

##########################################
## Plot cluster 12 SNP tree (Figure 3C) ##
##########################################

network12 <- ggtree(network12_tree)%<+% APHA_metadata_finallocations + 
  theme(legend.position="right") + 
  geom_tippoint(aes(colour = host), size=1.5) +
  scale_color_manual(values = c("#ff0000", "#0000cc")) +
  labs(color = "Host animal") +
  geom_tiplab(aes(subset = label %in% hubs, label = animalid), align = T, size =3, hjust = -0.2)+
  theme_tree2()+
  xlab("SNP distance") +
  theme(text = element_text(size = 15), axis.text.x = element_text(colour = "black")) +
  scale_x_continuous(breaks = c(0,2,4,6,8,10,12), labels = c(0,2,4,6,8,10,12)) +
  xlim_tree(14)

######################################
## Add trial area as metadata strip ##
######################################

network12Meta <- gheatmap(network12, network12All, width = 0.1, colnames = F, offset = 2, color = NULL) +
  scale_fill_manual(values = trial_id_colours, na.translate=FALSE, name = "Trial area")

##################################################
## Add life history for each animal to SNP tree ##
##################################################

network12Meta2 <- facet_plot(network12Meta, panel = 'Life history', data = inferred_infection_time_12,
                             geom = geom_segment, 
                             mapping = aes(x = LQ, xend = LQ, y = y-0.3, yend = y+0.3),size = 1, colour = "green") 

network12Meta3 <- facet_plot(network12Meta2, panel = 'Life history', data = inferred_infection_time_12,
                             geom = geom_segment, 
                             mapping = aes(x = UQ, xend = UQ, y = y-0.3, yend = y+0.3),size = 1, colour = "green")

network12Meta4 <- facet_plot(network12Meta3, panel = 'Life history', data = death_date_12,
                             geom = geom_segment, 
                             mapping = aes(x = deathdate, xend = deathdate, y = y-0.3, yend = y+0.3),size = 1, colour = "red")

network12Meta5 <- facet_plot(network12Meta4, panel = 'Life history', data = badgers_12,
                             geom = geom_segment, 
                             mapping = aes(x = badger_date, xend = badger_date, y = y-0.3, yend = y+0.3),size = 1, colour = "red")

network12Meta6 <- facet_plot(network12Meta5, panel = 'Life history', data = missing_animal_12,
                             geom = geom_segment, 
                             mapping = aes(x = deathdate, xend = deathdate-2, y = y, yend = y, colour = as.factor(onlocationkey)),
                             size = 1, linetype = "dotted", show.legend = FALSE) +
  scale_colour_manual(values = myColors, name = "Locations", limits = c("A","B", "C", "D", "E"))

network12Meta7 <- facet_plot(network12Meta6, panel = 'Life history', data = cattle_movements_12,
                             geom = geom_segment, 
                             mapping = aes(x = movementdate, xend = leaving_date, y = y, yend = y, colour = as.factor(onlocationkey)),
                             size = 1) +
  scale_colour_manual(values = myColors, name = "Locations", limits = c("A","B", "C", "D", "E")) +
  scale_x_continuous(breaks = c(1985,1990,1995,2000,2005,2010), labels = c(1985,1990,1995,2000,2005,2010)) +
  theme_bw()

#####################################
## Cluster 12 transmission network ##
#####################################

##################################
## Extract nodes for cluster 12 ##
##################################

Network12_nodes <- transNodes$network12 

#######################################
## Add hub labels to nodes dataframe ##
#######################################

Network12_nodes <- Network12_nodes %>% 
  mutate(hubs = ifelse(label %in% c("1", "2"), animalid,"")) 

Network12_nodes$onlocationkey <- Network12_nodes$onlocationkey %>% replace_na("No location")

#################################
## Create transmission network ##
#################################

Network12_TN <- graph_from_data_frame(transEdgesFilesRead$network12, vertices = Network12_nodes) %>% as_tbl_graph()

####################################
## Colour nodes based on location ##
####################################

myColors_TN <- c(myColors, "No location" = "grey") 

######################################################
## Plot cluster 12 transmission network (Figure 3D) ##
######################################################

plot_network12 <- ggraph(Network12_TN, layout = "nicely") +
  geom_edge_link(edge_colour = 'black', edge_alpha = 0.8, edge_width = 0.4,
                 arrow = arrow(length = unit(1, 'mm'), type = "closed"), end_cap = circle(1, "mm")) +
  geom_node_point(aes(shape = host, colour = as.factor(onlocationkey)), size = 4, show.legend= FALSE) + 
  scale_colour_manual(values = myColors_TN, name = "Location")+
  theme_graph()+
  geom_node_text(aes(label = hubs), repel = TRUE, size = 4, box.padding = 1)

################################
## Long distance transmission ##
################################

###############
## Cluster 6 ##
###############

#########################################
## Extract cluster 6 sequenced animals ##
#########################################

network6Samples <- all_sequenced_animals %>%
  left_join(complete_new_networks, by = c("LaneID" = "lane_id")) %>% 
  filter(new_network_id == "6") 

#################################################################
## Prepare dataframe for adding metadata to cluster 6 SNP tree ##
#################################################################

network6SamplesDF <- as.data.frame(network6Samples[,5])

rownames(network6SamplesDF) = network6Samples$LaneID

#########################################
## Plot cluster 6 SNP tree (Figure 3E) ##
#########################################

network6_tree_movement <- ggtree(network6_tree)%<+% APHA_metadata + 
  theme(legend.position="right") + 
  geom_tippoint(aes(color = host), size=1.5) +
  scale_color_manual(values = c("#ff0000", "#0000cc")) +
  geom_tiplab(aes(subset = (animalid %in% c(1, 2, 3, 4)), label = animalid), align = T, size =4,hjust = -0.1) +
  labs(color = "Host animal") +
  theme_tree2() +
  xlab("SNP distance") +
  theme(text = element_text(size = 15), axis.text.x = element_text(colour = "black")) + 
  xlim_tree(16)

######################################
## Add trial area as metadata strip ##
######################################

network6_tree_movementMeta <- gheatmap(network6_tree_movement, network6SamplesDF, width = 0.1, colnames = F, offset = 5, color = NULL) +
  scale_fill_manual(values = trial_id_colours, na.translate=FALSE, name = "Trial area")

##############################################################################
## Create map of cluster 6 isolates and add movements for selected isolates ##
##############################################################################

####################################################
## Extract movement data for animals in cluster 6 ##
####################################################

cattle_movements_network6 <- seq_animal_moves_continuous_7_days %>% 
  left_join(bovine_identities[,c("animalid", "new_network_id")], by = "animalid") %>% 
  filter(new_network_id == "6") %>% 
  filter(!is.na(onlocationkey) & !is.na(offlocationkey)) %>% 
  filter(onlocationkey != offlocationkey) %>% 
  arrange(animalid, desc(movementid)) %>% 
  left_join(location_Lat_Longs, by = c("offlocationkey" = "locationkey")) %>% 
  dplyr::rename(off_Long = Long, off_Lat = Lat, off_fuzzed_Long = fuzzed_Long, off_fuzzed_Lat = fuzzed_Lat) %>% 
  left_join(location_Lat_Longs, by = c("onlocationkey" = "locationkey")) %>% 
  dplyr::rename(on_Long = Long, on_Lat = Lat, on_fuzzed_Long = fuzzed_Long, on_fuzzed_Lat = fuzzed_Lat) %>% 
  add_count(offlocationkey, onlocationkey)

#################################################
## Specify geographic coordinates for map plot ##
#################################################

Long_max <- -4.2
Long_min <- -4.75
Lat_min <- 50.3
Lat_max <- 50.9

################################################
## Extract movement data for selected animals ##
################################################

plot6_moves <- cattle_movements_network6 %>% 
  filter(animalid %in% c(1, 2, 3, 4))

#####################################################
## Extract final locations for cattle in cluster 6 ##
#####################################################

bovine_locations_6 <- animal_final_locations %>% 
  filter(host == "Bovine") %>% 
  filter(new_network_id =="6") %>% 
  add_count(onlocationkey)

######################################################
## Extract final locations for badgers in cluster 6 ##
######################################################

badger_locations6 <- animal_final_locations %>%
  filter(host == "Badger") %>% 
  filter(new_network_id == "6") 

########################################
## Create labels for cattle locations ##
########################################

bov_loc_labs <- bovine_locations_6 %>% 
  select(onlocationkey, Long, Lat, fuzzed_Long, fuzzed_Lat) %>% 
  unique()

##################################################
## Extract location labels for selected animals ##
##################################################

plot6_bovlocs <- bov_loc_labs %>% 
  filter(onlocationkey %in% c("A", "B", "C", "D"))

######################################################
## Extract location coordinates for missing animals ##
######################################################

missing_animals_6_latlong <- fuzzed_missing_LatLong %>% 
  filter(new_network_id == "6") %>% 
  add_count(locationkey)

################################################
## Plot map of cluster 6 isolates (Figure 3F) ##
################################################

plot6_geo <- ggplot() + 
  geom_path(data = UK, aes(x = long, y = lat, group = group)) +
  coord_map(xlim = c(Long_min, Long_max), ylim = c(Lat_min, Lat_max)) +
  xlim(c(Long_min, Long_max))+ # remove to include those going outside
  ylim(c(Lat_min, Lat_max))+ # remove to include those going outside
  theme_bw()+
  geom_point(data = bovine_locations_6, aes(x = fuzzed_Long, y = fuzzed_Lat, shape = host, size = n, colour = as.factor(onlocationkey)), stroke = 2, alpha = 0.5, show.legend = FALSE)+
  geom_point(data = missing_animals_6_latlong, aes(x = fuzzed_Long, y = fuzzed_Lat, shape = host, size = n, colour = as.factor(locationkey)), show.legend = FALSE)+
  scale_colour_manual(values = colourList, name = "Locations")+
  new_scale_colour()+
  geom_point(data = badger_locations6, aes(x = fuzzed_Long, y = fuzzed_Lat, shape = host), stroke = 1)+
  scale_shape_manual(values = shapes_animals)+
  labs(size = "No. of animals") +
  new_scale_colour()+
  geom_segment(data = plot6_moves, aes(x = off_fuzzed_Long, y = off_fuzzed_Lat, xend = on_fuzzed_Long, yend = on_fuzzed_Lat, colour = as.factor(animalid)), 
               arrow = arrow(length = unit(2, 'mm'), type = "open"))+
  labs(colour = "Bovine")+
  xlab("Longitude") + 
  ylab("Latitude")+
  geom_text_repel(data = plot6_bovlocs, aes(x = fuzzed_Long, y = fuzzed_Lat, label = onlocationkey), size = 4)+
  theme(axis.title=element_text(size=15)) + 
  theme(axis.text.x = element_text(size = 15)) +
  theme(axis.text.y = element_text(size = 15)) +
  ggsn::scalebar(UK, dist = 10, dist_unit = "km", transform = TRUE, model = "WGS84", 
                 st.size = 4, st.dist = 0.0015 , height = 0.001, anchor = c(x = Long_max -0.05, y = Lat_min + 0.03), location = "bottomright")
