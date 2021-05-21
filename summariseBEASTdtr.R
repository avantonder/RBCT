library(tidyverse)

#############################
### Set working directory ###
#############################

DTRDir <- "~/BEAST_DTR/"
setwd(DTRDir)

###################################
## List infection time csv files ##
###################################

DTRFiles <- list.files(pattern = "*stats.csv")

###################################
## Extract network names as list ##
###################################

DTRFilesNames <- gsub('_clock.rate.stats.csv', '', DTRFiles)

DTRFilesNames <- gsub('network', '', DTRFilesNames)

#######################
## Read in csv files ##
#######################

DTRFilesRead <- lapply(DTRFiles, function(x){
  read_csv(x)
})

names(DTRFilesRead) <- DTRFilesNames

###############################################
## Collapse list of DTR stats into single df ##
###############################################

DTRdf <- bind_rows(DTRFilesRead, .id = "network_id")

######################################################################################################
## Create plot of clock rates for real and randomized runs for each network (Supplementary Figure 3 ##
######################################################################################################

DTRdfColour <- DTRdf %>%
  group_by(network_id) %>%
  mutate(colour = ifelse(calibr == 0, paste("red"), paste("black")))

DTRdfLines <- DTRdfColour[DTRdfColour$calibr == 0,]

DTRplot <- ggplot(data = transform(DTRdfColour,
                                   network_id = factor(network_id,
                                                       levels = c(1,2,3,4,5,6,7,8,9,10,11,12)))) + 
  geom_point(aes(x = factor(calibr), y = median, color = colour)) + 
  geom_errorbar(aes(x = factor(calibr), ymin = lowerHPD, ymax = HigherHPD, color = colour), 
                width = 0.2, position = position_dodge(0.9)) + 
  theme_bw() + 
  ylab("Clock rate") + 
  xlab("") +
  geom_hline(data = transform(DTRdfLines,
                              network_id = factor(network_id,
                                                  levels = c(1,2,3,4,5,6,7,8,9,10,11,12))), 
             aes(yintercept = lowerHPD), colour = "red", linetype = "dashed") +
  geom_hline(data = transform(DTRdfLines,
                              network_id = factor(network_id,
                                                  levels = c(1,2,3,4,5,6,7,8,9,10,11,12))), 
             aes(yintercept = HigherHPD), colour = "red", linetype = "dashed") +
  theme(axis.text.x = element_blank()) +
  facet_wrap(~network_id, ncol = 4) +
  scale_color_manual(values = c("black", "red")) + 
  theme(strip.text.x = element_text(size = 15, color = "black"),
        strip.background = element_rect(color="black", fill="white", size=1, linetype="solid")) +
  theme(axis.title=element_text(size=15)) + 
  theme(axis.text.y = element_text(size = 15)) +
  theme(legend.position = "none")