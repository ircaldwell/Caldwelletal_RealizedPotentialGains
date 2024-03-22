###### Caldwell et al. Code for "Protection efforts have resulted in ~10% of existing fish biomass on global coral reefs" ####
######  Create main figures for paper ####
######  Code author: Iain R. Caldwell
######  Date last revised: March 21, 2024
######  Main figures for paper:
######    Fig 1. Realized gains
######      a. Cumulative biomass vs. % of sites for status quo and fished scenario
######      b. Realized gains given subsampled full protection
######      c. Map of realized gains
######      Assemble figure 1
######    Fig 2. Variaion in realized gains within protection categories
######      a. Stacked distribution of realized gains (colored by management)
######      b. Realized coral reef fish biomass gains along a gradient of nearest market gravity
######      Assemble figure 2
######    Fig 3. Potential gains
######      a. Cumulative potential (positive) and realized (negative) gains in order but colored by management 
######      e. Map of potential full protection gains  
######      Assemble figure 3
######  Supplementary figures for paper:
######    Fig S1. Coral reef fish biomass predicted by socio-ecological context
######      a. Map of surveyed sites sized by status quo biomass and colored by protection category
######      b. Stacked histogram of status quo fish biomass
######      c. Effect sizes for all fixed covariates in the most predictive spatial GLMMs
######    Fig S2. Comparison of residuals with distance to the nearest MPA ####
######    Fig S3. Test of size and age breaks to differentiate big and old vs. small or new fully protected MPAs
######    Fig S4. Random effects in the most predictive spatial GLMMs
######    Fig S5. Model performance/diagnostics
######    Fig S6. Representativeness of the most predictive social and environmental predictors, ordered by absolute effect size 
######    Fig S7. Representativeness of social and environmental predictors not in one of the best models

rm(list = ls()) #remove past stored objects
options(scipen=999) #disable scientific notation
set.seed(1234)

#Load libraries
library(ggplot2) #plotting
library(tidyverse) 
library(ggpubr) #combining plots
library(lemon) #for facet_rep_grid function
library(plotrix) #for calculating uncertainty (SE) among sites 
library(ggridges) #ridge plots
library(scales) #for removing trailing zeroes in log plots and custom transformations
library(spaMM) #for running the spatial mixed models
library(mgcv) #for gam in Fig S2
library(gratia) #for testing gams in Fig S2
library(emmeans) #estimating marginal means
library(DHARMa) #for model diagnostics
library(gap) #for model diagnostics in DHARMa

#Set random seed for reproducibility
set.seed(1234)

#Set the number of samples
numSamples = 1000

####Load functions ####
#To darken colors
darken <- function(color, factor=1.4){
  col <- col2rgb(color)
  col <- col/factor
  col <- rgb(t(col), maxColorValue=255)
  col
}

####Set parameters and directories ####
dataDir <- "../Data/"
modelDir <- "../Models/"
plotDir <- "../Plots/"

#Check if the plot folder exists and create it if not
if(!file.exists(plotDir)) {
  dir.create(file.path(plotDir))
} 

######  Load predictions from "Caldwelletal_RealizedPotentialGains_PredictGains.R" code ####
###Open RDS file with predictions, best model, and reef mask data (for supplementary figure)
predictFilename <- paste0(dataDir, "Caldwelletal_RealizedPotentialGains_FishBiomassGainsPredSimsCovariates.rds")
if(file.exists(predictFilename)) {
  predictTBL <- readRDS(file = predictFilename)
} else {
  message("Predict file missing - running predict code")
  source("Caldwelletal_RealizedPotentialGains_PredictGains.R")
}

#Get numbers for each protection type
manageSamples <- as.data.frame(table(predictTBL$Management))

######  Load list with best models ####
bestSpammResTBL <- read_csv(file = paste0(dataDir, "Caldwelletal_RealizedPotentialGains_BestSpammSummResWithFilenames.csv")) 

#get the number of best models
numModels <- nrow(bestSpammResTBL)

######  Load the reef mask data (for supplement) ####
reefMaskTBL <- read_csv(file = paste0(dataDir, "Caldwelletal_RealizedPotentialGains_ReefMaskCovariates.csv"))

######    Fig 1. Realized gains ####
######      > 1a. Cumulative biomass vs. % of sites for status quo and fished scenario ####
###Set colors for management types
manageColors <- c("#03018C", "#4259C3", "#9EC2FF", "#009E73", "#E69F00")
names(manageColors) <- c("UnfishedHighBigOld", "UnfishedHighSmallNew", "UnfishedLow", "Restricted", "Fished")

###Set labels for management types
manageLabels <- c("Big and old high compliance\nfully protected MPAs",
                  "Small or new high compliance\nfully protected MPAs",
                  "Low compliance\nfully protected MPAs",
                  "Restricted fishing",
                  "Open to fishing")
names(manageLabels) <- c("UnfishedHighBigOld", "UnfishedHighSmallNew", "UnfishedLow", "Restricted", "Fished")

###Set new labels with sample sizes
manageLabelsWithSampleSizes <- c(paste0("Big and old high compliance\nfully protected MPAs (", manageSamples$Freq[manageSamples$Var1 == "UnfishedHighBigOld"], " sites)"),
                                 paste0("Small or new high compliance\nfully protected MPAs (", manageSamples$Freq[manageSamples$Var1 == "UnfishedHighSmallNew"], " sites)"),
                                 paste0("Low compliance\nfully protected MPAs (", manageSamples$Freq[manageSamples$Var1 == "UnfishedLow"], " sites)"),
                                 paste0("Restricted fishing\n(", manageSamples$Freq[manageSamples$Var1 == "Restricted"], " sites)"),
                                 paste0("Open to fishing\n(", manageSamples$Freq[manageSamples$Var1 == "Fished"], " sites)"))
names(manageLabelsWithSampleSizes) <- c("UnfishedHighBigOld", "UnfishedHighSmallNew", "UnfishedLow", "Restricted", "Fished")

###Set colors for gains
gainColors <- c("#dfe0e2", "#5ba56e", "#236477")
names(gainColors) <- c("RealizedGains", "RestrictedGains", "PotentialGains")

###Sort data for figure and calculate the cumulative biomasses as percentages of status quo
cumBiomassTBL <- predictTBL %>% 
  select(Management, contains("StatusQuoPercTotalStatusQuoSample_"), contains("FishedPercTotalStatusQuoSample_"), StatusQuoBiomassMed) %>% 
  mutate(Management = factor(Management, levels = c("UnfishedHighBigOld", "UnfishedHighSmallNew", "UnfishedLow", "Restricted", "Fished"))) %>% 
  arrange(desc(Management), desc(StatusQuoBiomassMed)) %>% 
  mutate(PercentageSites = c(1:nrow(predictTBL))/nrow(predictTBL)*100)
    
cumStatusQuoPercStatusQuoQuantiles <- select(cumBiomassTBL, starts_with("StatusQuoPercTotalStatusQuoSample_")) %>%
  mutate(across(everything(), ~cumsum(.))) %>%
  rename_with(~paste0(., "_cumulative"), everything()) %>%
  rowwise() %>%
  mutate(
    CumStatusQuoPercTotalStatusQuoMed = median(c_across(starts_with("StatusQuoPercTotalStatusQuoSample_")), na.rm = TRUE),
    CumStatusQuoPercTotalStatusQuoQuant_025 = quantile(c_across(starts_with("StatusQuoPercTotalStatusQuoSample_")), probs = 0.025, na.rm = TRUE),
    CumStatusQuoPercTotalStatusQuoQuant_975 = quantile(c_across(starts_with("StatusQuoPercTotalStatusQuoSample_")), probs = 0.975, na.rm = TRUE)
  ) %>% 
  select(CumStatusQuoPercTotalStatusQuoMed, CumStatusQuoPercTotalStatusQuoQuant_025, CumStatusQuoPercTotalStatusQuoQuant_975)

cumFishedPercStatusQuoQuantiles <- select(cumBiomassTBL, starts_with("FishedPercTotalStatusQuoSample_")) %>%
  mutate(across(everything(), ~cumsum(.))) %>%
  rename_with(~paste0(., "_cumulative"), everything()) %>%
  rowwise() %>%
  mutate(
    CumFishedPercTotalStatusQuoMed = median(c_across(starts_with("FishedPercTotalStatusQuoSample_")), na.rm = TRUE),
    CumFishedPercTotalStatusQuoQuant_025 = quantile(c_across(starts_with("FishedPercTotalStatusQuoSample_")), probs = 0.025, na.rm = TRUE),
    CumFishedPercTotalStatusQuoQuant_975 = quantile(c_across(starts_with("FishedPercTotalStatusQuoSample_")), probs = 0.975, na.rm = TRUE)
  ) %>% 
  select(CumFishedPercTotalStatusQuoMed, CumFishedPercTotalStatusQuoQuant_025, CumFishedPercTotalStatusQuoQuant_975)

cumBiomassTBL <- bind_cols(cumBiomassTBL, cumStatusQuoPercStatusQuoQuantiles, cumFishedPercStatusQuoQuantiles)

#### Get numbers for paper ####
#Percentage realized gains from high compliance MPAs
realizedGainsHighCompMPaPercStatusQuoSums <- predictTBL %>% 
  filter(Management %in% c("UnfishedHighSmallNew", "UnfishedHighBigOld")) %>% 
  select(starts_with("RealizedGainsPercStatusQuoSample_")) %>%
  summarize(across(everything(), ~sum(.)))

median(as.numeric(realizedGainsHighCompMPaPercStatusQuoSums)) #14.2% --> realized gains from high compliance MPAs
quantile(as.numeric(realizedGainsHighCompMPaPercStatusQuoSums), probs = 0.025) #13.5
quantile(as.numeric(realizedGainsHighCompMPaPercStatusQuoSums), probs = 0.975) #15.0

realizedGainsAllPercStatusQuoSums <- predictTBL %>% 
  select(starts_with("RealizedGainsPercStatusQuoSample_")) %>%
  summarize(across(everything(), ~sum(.)))
  
realizedGainsPercHighCompMpaVsAllMedQuant <- as.numeric(realizedGainsHighCompMPaPercStatusQuoSums/realizedGainsAllPercStatusQuoSums)*100 

median(realizedGainsPercHighCompMpaVsAllMedQuant) #66.6
quantile(realizedGainsPercHighCompMpaVsAllMedQuant, probs = 0.025) #64.1
quantile(realizedGainsPercHighCompMpaVsAllMedQuant, probs = 0.975) #69.0

realizedGainsRestPercStatusQuoSums <- predictTBL %>% 
  filter(Management %in% c("Restricted")) %>% 
  select(starts_with("RealizedGainsPercStatusQuoSample_")) %>%
  summarize(across(everything(), ~sum(.)))

median(as.numeric(realizedGainsRestPercStatusQuoSums)) #5.7% --> realized gains from high compliance MPAs
quantile(as.numeric(realizedGainsRestPercStatusQuoSums), probs = 0.025) #5.2
quantile(as.numeric(realizedGainsRestPercStatusQuoSums), probs = 0.975) #6.3

###>> Fig 1a plot ####
Fig1a_CumBiomassVsPerSitesFishedStatusQuoRealGainsPlot <- ggplot(data = cumBiomassTBL, aes(x = PercentageSites)) +
  geom_ribbon(aes(ymin = max(CumFishedPercTotalStatusQuoQuant_025), ymax = max(CumFishedPercTotalStatusQuoQuant_975)),
              fill = gainColors["RealizedGains"]) +
  geom_segment(
    aes(y = max(CumFishedPercTotalStatusQuoMed), yend = max(CumFishedPercTotalStatusQuoMed)), x = 0, xend = 100, 
    linewidth = 1, 
    linetype = "dashed",
    colour = "black" 
  ) +
  geom_segment(
    y = 100, yend = 100, x = 0, xend = 100, 
    linewidth = 1, 
    linetype = "dashed",
    colour = "black" 
  ) +
  geom_ribbon(data = cumBiomassTBL %>% filter(Management != "Fished"),
              aes(ymin = CumFishedPercTotalStatusQuoQuant_025,
                  ymax = CumFishedPercTotalStatusQuoQuant_975),
              fill = manageColors["Fished"], alpha = 0.5) +
  geom_ribbon(aes(ymin = CumStatusQuoPercTotalStatusQuoQuant_025,
                  ymax = CumStatusQuoPercTotalStatusQuoQuant_975,
                  fill = Management), alpha = 0.5) +
  geom_line(data = cumBiomassTBL %>% filter(Management != "Fished"),
            aes(y = CumFishedPercTotalStatusQuoMed),
            color = manageColors["Fished"], linewidth = 1, linetype = "dashed") +
  geom_line(aes(y = CumStatusQuoPercTotalStatusQuoMed, color = Management), linewidth = 1) +
  scale_x_continuous(limits = c(0,101),
                     expand = c(0, 0),
                     name = "% of sites",
                     breaks = c(0, 20, 40, 60, 80, 100)) +
  scale_y_continuous(limits = c(0, 101),
                     expand = c(0,0),
                     name = "Cumulative fish biomass\n(% of status quo total)",
                     breaks = c(0, 20, 40, 60, 80, 100)) +
  scale_color_manual(values = manageColors,
                     labels = manageLabels,
                     #name = "Current (status quo) management:",
                     name = NULL) +
  scale_fill_manual(values = manageColors,
                     labels = manageLabels,
                     #name = "Current (status quo) management:",
                     name = NULL) +
  theme_bw() +
  guides(color = "none", fill = "none") +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = 'black'),
        aspect.ratio = 1) +
  #Add label for realized gains
  annotate(geom = "text", x = 79, y = 90, label = "total\nrealized\ngains", color = "black",
           hjust = 1, vjust = 0.5) +
  #Add arrow for total realized gains
  geom_segment(
    aes(x = 80, y = max(CumFishedPercTotalStatusQuoMed)), xend = 80, yend = 100,
    lineend = "round", 
    linejoin = "round",
    linewidth = 1, 
    arrow = arrow(length = unit(0.3, "cm")),
    colour = "black" 
  ) +
  #Add arrow and text for the status quo biomass
  geom_segment(
    aes(x = 59, y = 67, xend = 64, yend = 62),
    lineend = "round", 
    linejoin = "round",
    linewidth = 1, 
    arrow = arrow(length = unit(0.3, "cm")),
    colour = "black" 
  ) + 
  annotate(geom = "text", x = 50, y = 70, label = "Status quo", color = "black") +
  #Add arrow pointing to the fished scenario
  geom_segment(
    aes(x = 80, y = 50, xend = 75, yend = 59),
    lineend = "round", 
    linejoin = "round",
    linewidth = 1, 
    arrow = arrow(length = unit(0.3, "cm")),
    colour = "black" 
  ) + 
  annotate(geom = "text", x = 85, y = 48, label = "Fished scenario", color = "black")

######      > 1b. Realized gains given subsampled full protection ####
# Create new dataset for figure with actual data and subsampled data (1000 subsamples) --> 2.8% would be 40 (Full), leaving 74 Partial
#Get the percentages for the actual data
realGainsFullMpaByPercFilename <- paste0(dataDir, "Caldwelletal_RealizedPotentialGains_SubsampledByPercFullMPA.csv")

if(file.exists(realGainsFullMpaByPercFilename)) {
  message("Realized gains already subsampled")
  realGainsFullMpaByPercProtTBL <- read_csv(realGainsFullMpaByPercFilename)
} else {
  realGainsFullMpaPercStatusQuoActual <- c()
  
  i = 1
  for(i in 1:numSamples) {
    realGainsFullMpaPercStatusQuoActual <- c(realGainsFullMpaPercStatusQuoActual,
                                             sum(predictTBL[,paste0("RealizedGainsSample_",i)])/
                                               sum(predictTBL[,paste0("StatusQuoPredSample_",i)]))
  }
  
  realGainsFullMpaByPercProtTBL <- data.frame(NumFullMPA = sum(!predictTBL$Management %in% c("Fished", "Restricted")),
                                              PercFullMPA = sum(!predictTBL$Management %in% c("Fished", "Restricted"))/nrow(predictTBL)*100,
                                              RealGainsFullMpaPercStatusQuoMed = unname(median(x = realGainsFullMpaPercStatusQuoActual)),
                                              RealGainsFullMpaPercStatusQuoQuant975 = unname(quantile(x = realGainsFullMpaPercStatusQuoActual, probs = 0.975)),
                                              RealGainsFullMpaPercStatusQuoQuant025 = unname(quantile(x = realGainsFullMpaPercStatusQuoActual, probs = 0.025)),
                                              Data = "All surveys")
  
  fullProtTBL <- predictTBL %>% 
    filter(!Management %in% c("Fished", "Restricted"))
  
  noFullProtTBL <- predictTBL %>% 
    filter(Management %in% c("Fished", "Restricted"))
  
  #Add another row to the dataframe for no MPAs
  realGainsFullMpaByPercProtTBL[nrow(realGainsFullMpaByPercProtTBL)+1,] <- NA
  realGainsFullMpaByPercProtTBL$NumFullMPA[nrow(realGainsFullMpaByPercProtTBL)] <- 0
  realGainsFullMpaByPercProtTBL$PercFullMPA[nrow(realGainsFullMpaByPercProtTBL)] <- 0
  
  realGainsFullMpaPercStatusQuoGlobalProt <- c()
  
  j = 1
  for(j in 1:numSamples) {
    realGainsFullMpaPercStatusQuoGlobalProt <- c(realGainsFullMpaPercStatusQuoGlobalProt,
                                                 sum(noFullProtTBL[,paste0("RealizedGainsSample_",j)])/
                                                   sum(noFullProtTBL[,paste0("StatusQuoPredSample_",j)]))
  }
  
  realGainsFullMpaByPercProtTBL$RealGainsFullMpaPercStatusQuoMed[nrow(realGainsFullMpaByPercProtTBL)] <- unname(median(x = realGainsFullMpaPercStatusQuoGlobalProt))
  realGainsFullMpaByPercProtTBL$RealGainsFullMpaPercStatusQuoQuant975[nrow(realGainsFullMpaByPercProtTBL)] <- unname(quantile(x = realGainsFullMpaPercStatusQuoGlobalProt, probs = 0.975))
  realGainsFullMpaByPercProtTBL$RealGainsFullMpaPercStatusQuoQuant025[nrow(realGainsFullMpaByPercProtTBL)] <- unname(quantile(x = realGainsFullMpaPercStatusQuoGlobalProt, probs = 0.025))
  realGainsFullMpaByPercProtTBL$Data[nrow(realGainsFullMpaByPercProtTBL)] <- "No Full MPAs"
  
  #Successively subsample the MPAs and recalculate the realized gains
  subsamples <- seq(from = 1, to = sum(predictTBL$Management %in% c("UnfishedHighBigOld", "UnfishedHighSmallNew", "UnfishedLow"))-1, by = 5)
  
  i = 1
  for(i in subsamples) {
    message("Started subsampling ", i, " full MPAs")
    #Check if this number of full MPAs has already been run
    if(i %in% realGainsFullMpaByPercProtTBL$NumFullMPA) {
      message("This subsample has already been run")
    } else {
      realGainsFullMpaByPercProtTBL[nrow(realGainsFullMpaByPercProtTBL)+1,] <- NA
      realGainsFullMpaByPercProtTBL$NumFullMPA[nrow(realGainsFullMpaByPercProtTBL)] <- i
      realGainsFullMpaByPercProtTBL$PercFullMPA[nrow(realGainsFullMpaByPercProtTBL)] <- i/(nrow(noFullProtTBL)+i)*100
      
      realGainsFullMpaPercStatusQuoGlobalProt <- c()
      
      j = 1
      for(j in 1:numSamples) {
        fullRowSample <- sample(x = row.names(fullProtTBL), size = i)
        subsampleTBL <- bind_rows(fullProtTBL[fullRowSample,],
                                  noFullProtTBL)
        realGainsFullMpaPercStatusQuoGlobalProt <- c(realGainsFullMpaPercStatusQuoGlobalProt,
                                                     sum(subsampleTBL[,paste0("RealizedGainsSample_",j)])/
                                                       sum(subsampleTBL[,paste0("StatusQuoPredSample_",j)]))
      }
      
      
      realGainsFullMpaByPercProtTBL$RealGainsFullMpaPercStatusQuoMed[nrow(realGainsFullMpaByPercProtTBL)] <- unname(median(x = realGainsFullMpaPercStatusQuoGlobalProt))
      realGainsFullMpaByPercProtTBL$RealGainsFullMpaPercStatusQuoQuant975[nrow(realGainsFullMpaByPercProtTBL)] <- unname(quantile(x = realGainsFullMpaPercStatusQuoGlobalProt, probs = 0.975))
      realGainsFullMpaByPercProtTBL$RealGainsFullMpaPercStatusQuoQuant025[nrow(realGainsFullMpaByPercProtTBL)] <- unname(quantile(x = realGainsFullMpaPercStatusQuoGlobalProt, probs = 0.025))
      realGainsFullMpaByPercProtTBL$Data[nrow(realGainsFullMpaByPercProtTBL)] <- "Subsample"
    }
  }
  
  realGainsFullMpaByPercProtTBL <- realGainsFullMpaByPercProtTBL %>% 
    filter(!is.na(RealGainsFullMpaPercStatusQuoMed))
  
  #Save the dataset so this doesn't have to be run again
  write_csv(x = realGainsFullMpaByPercProtTBL, file = realGainsFullMpaByPercFilename)
}

realGainsFullMpaByPercProtTBL <- realGainsFullMpaByPercProtTBL %>% 
  arrange(NumFullMPA)

#Get the closest percentage to the Allen Coral percentages (3.06%)
closestPercAllenCoralProt <- realGainsFullMpaByPercProtTBL$PercFullMPA[which.min(abs(realGainsFullMpaByPercProtTBL$PercFullMPA - 3.06))] 
realGainsFullMpaByPercProtTBL$NumFullMPA[which.min(abs(realGainsFullMpaByPercProtTBL$PercFullMPA - 3.06))]#61

#Get the median and quantiles for realized gains if we had sampled 3% MPAs
realGainsFullMpaByPercProtTBL %>% 
  filter(NumFullMPA == 61) #10.2% (9.4 - 11.1%)

###>> Fig 1b plot ####
Fig1b_RealGainsSubsampleAllenCoralPercPlot <- ggplot(data = realGainsFullMpaByPercProtTBL,
                                               aes(x = PercFullMPA,
                                                   y = RealGainsFullMpaPercStatusQuoMed*100,
                                                   ymin = RealGainsFullMpaPercStatusQuoQuant025*100,
                                                   ymax = RealGainsFullMpaPercStatusQuoQuant975*100)) +
  geom_line(size = 1) +
  geom_ribbon(alpha = 0.5) +
  geom_segment(data = realGainsFullMpaByPercProtTBL %>% filter(Data == "All surveys"),
               aes(x = 0, xend = PercFullMPA, y = RealGainsFullMpaPercStatusQuoMed*100, yend = RealGainsFullMpaPercStatusQuoMed*100),
               linetype = "dashed") +
  geom_segment(data = realGainsFullMpaByPercProtTBL %>% filter(Data == "All surveys"),
               aes(x = PercFullMPA, xend = PercFullMPA, y = 5, yend = RealGainsFullMpaPercStatusQuoMed*100),
               linetype = "dashed") +
  geom_errorbar(data = realGainsFullMpaByPercProtTBL %>% filter(Data == "All surveys")) +
  geom_point(data = realGainsFullMpaByPercProtTBL %>% filter(Data == "All surveys"), size = 3, shape = 21, fill = "grey") +
  geom_segment(data = realGainsFullMpaByPercProtTBL %>% filter(PercFullMPA == closestPercAllenCoralProt),
               aes(x = 0, xend = PercFullMPA, y = RealGainsFullMpaPercStatusQuoMed*100, yend = RealGainsFullMpaPercStatusQuoMed*100),
               linetype = "dotted") +
  geom_segment(data = realGainsFullMpaByPercProtTBL %>% filter(PercFullMPA == closestPercAllenCoralProt),
               aes(x = PercFullMPA, xend = PercFullMPA, y = 5, yend = RealGainsFullMpaPercStatusQuoMed*100),
               linetype = "dotted") +
  geom_errorbar(data = realGainsFullMpaByPercProtTBL %>% filter(PercFullMPA == closestPercAllenCoralProt)) +
  geom_point(data = realGainsFullMpaByPercProtTBL %>% filter(PercFullMPA == closestPercAllenCoralProt), size = 3, shape = 21, fill = "white") +
  #geom_segment(x = 0, xend = 0, y = 0, yend = )
  scale_y_continuous(name = "Realized gains\n(% of total status quo biomass)",
                     #limits = c(5,20),
                     expand = c(0,0)
  ) +
  scale_x_continuous(name = "% fully protected MPA coverage",
                     #limits = c(0,26),
                     expand = c(0,0)) +
  theme_classic() +
  theme(aspect.ratio = 1) 

######      > 1c. Map of realized gains ####
###Set the minimum value for the log gains
zeroLogTrans = 1

###Set up map
mapWorld <- map_data('world2')

###Prepare data
realGainMapTBL <- predictTBL %>%
  select(Management, Long_DD, Lat_DD, RealizedGainMed) %>% 
  mutate(Long_DD = ifelse(test = Long_DD < -25, yes = Long_DD + 360, no = Long_DD),
         RealizedGainsLogPos = ifelse(test = RealizedGainMed <= zeroLogTrans, yes = zeroLogTrans, no = RealizedGainMed),
         Management = factor(Management, levels = c("UnfishedHighBigOld", "UnfishedHighSmallNew", "UnfishedLow", "Restricted", "Fished"))) %>% 
  arrange(desc(RealizedGainsLogPos))

#Create new gain limits
mapGain_breaks <- c(1, 10, 100, 1000)
mapGain_labels <- c(0, 10, 100, 1000) 
gainAmount_limits <- c(zeroLogTrans, max(predictTBL$RealizedGainHigh975))

####>> Fig 1c map ####
Fig1c_RealGainsMap_SizeLegendBelow <- ggplot() + 
  geom_polygon(data = mapWorld, aes(x = long, y = lat, group = group), fill = "grey50", color = "grey50") +
  coord_fixed(ratio = 1.3, xlim = c(40, 320), ylim = c(-25, 25)) +
  geom_jitter(data = realGainMapTBL, 
              aes(x = Long_DD, y = Lat_DD, fill = Management, size = RealizedGainsLogPos),
              colour = "black", pch = 21, width = 1, height = 1) +
  scale_fill_manual(name = NULL,
                    values = manageColors,
                    labels = manageLabelsWithSampleSizes) +
  scale_size_continuous(breaks = mapGain_breaks, 
                        labels = mapGain_labels, 
                        range = c(1,8), 
                        limits = gainAmount_limits) +
  guides(fill = "none",
         size = guide_legend(title = "Realized fish biomass gains (kg/ha):"), order = 2) +
  geom_hline(yintercept = 23.43695, lty = 2) +
  geom_hline(yintercept = -23.43695, lty = 2) + 
  scale_x_continuous("", breaks = c(30, 80, 130, 180, 230, 280, 330), labels=c(-150, -100, -50, 0, 50, 100, 150)) +
  scale_y_continuous("", breaks = c(-20, -10, 0, 10, 20)) +
  theme(plot.subtitle = element_text(face = "italic"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        #legend.box.spacing = unit(1, 'cm'),
        legend.position = "bottom",
        legend.key = element_rect(fill = NA),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-20, -20, -20, -20),
        legend.text.align = 0.5)

RealGainsMap_WithOnlyColorLegend <- ggplot() + 
  geom_polygon(data = mapWorld, aes(x = long, y = lat, group = group), fill = "grey50", color = "grey50") +
  coord_fixed(ratio = 1.3, xlim = c(40, 320), ylim = c(-25, 25)) +
  geom_jitter(data = realGainMapTBL, 
              aes(x = Long_DD, y = Lat_DD, fill = Management, size = RealizedGainsLogPos),
              colour = "black", pch = 21, width = 1, height = 1) +
  scale_fill_manual(name = NULL,
                    values = manageColors,
                    labels = manageLabels) +
  guides(fill = guide_legend(override.aes = list(shape = 22, size = 6), order = 1, reverse = T),
         size = "none") +
  theme(legend.box.spacing = unit(1, 'cm'),
        legend.position = "top",
        legend.box = "vertical",
        legend.key = element_rect(fill = NA),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-20, -20, -20, -20),
        legend.text.align = 0.5)

ProtectionCategory_FillLegend <- get_legend(RealGainsMap_WithOnlyColorLegend)

######      > Assemble figure 1 ####
Fig1_RealGainsCumBiomassCorrSubsampleMap <- ggpubr::ggarrange(
  #First row with cumulative realized gains and subsampled realized gains vs. fully protected MPA coverage 
  ggpubr::ggarrange(Fig1a_CumBiomassVsPerSitesFishedStatusQuoRealGainsPlot,
                    Fig1b_RealGainsSubsampleAllenCoralPercPlot,
                    ncol = 2, labels = c("a", "b"), align = "v"),
  #Second row with fill legend
  ProtectionCategory_FillLegend,
  ggpubr::ggarrange(Fig1c_RealGainsMap_SizeLegendBelow, labels = c("c")),
  nrow = 3, heights = c(1.5,0.1,1)
)

#Plot as a tiff
ggsave(plot = Fig1_RealGainsCumBiomassCorrSubsampleMap,
       filename = paste0(plotDir, "Caldwelletal_RealizedPotentialGains_Fig1_RealGainsCumBiomassCorrSampleMap.tiff"),
       dpi = 1000, device = "tiff", bg = "white", width = 10, height = 10)

#Plot as a PDF
ggsave(plot = Fig1_RealGainsCumBiomassCorrSubsampleMap,
       filename = paste0(plotDir, "Caldwelletal_RealizedPotentialGains_Fig1_RealGainsCumBiomassCorrSampleMap.pdf"),
       dpi = 1000, device = "pdf", bg = "white", width = 10, height = 10)


######    Fig 2. Variation in realized gains within protection categories ####
######      > 2a. Stacked distribution of realized gains (colored by management) ####
###Create a new tibble for the plot that changes 0 values to 1 for log plotting
stackRealTBL <- predictTBL %>% 
  mutate(RealizedGainsLogPos = ifelse(test = RealizedGainMed <= zeroLogTrans, yes = zeroLogTrans, no = RealizedGainMed),
         Management = factor(Management, levels = c("UnfishedHighBigOld", "UnfishedHighSmallNew", "UnfishedLow", "Restricted", "Fished")))

###Set the gain limits
gainAmount_breaks <- c(1, 10, 100, 1000)
gainAmount_labels <- c(0, 10, 100, 1000)

#Make a dataset with the medians and quantiles for each management category (from the samples)
restRG <- c()
unfishedLowRG <- c()
unfishedHighSmallNewRG <- c()
unfishedHighBigOldRG <- c()
allRG <- c()

for(i in 1:numSamples) {
  restRG <- c(restRG, unname(unlist(predictTBL[predictTBL$Management == "Restricted", paste0("RealizedGainsSample_", i)])))
  unfishedLowRG <- c(unfishedLowRG, unname(unlist(predictTBL[predictTBL$Management == "UnfishedLow", paste0("RealizedGainsSample_", i)])))
  unfishedHighSmallNewRG <- c(unfishedHighSmallNewRG, unname(unlist(predictTBL[predictTBL$Management == "UnfishedHighSmallNew", paste0("RealizedGainsSample_", i)])))
  unfishedHighBigOldRG <- c(unfishedHighBigOldRG, unname(unlist(predictTBL[predictTBL$Management == "UnfishedHighBigOld", paste0("RealizedGainsSample_", i)])))
  allRG <- c(allRG, unname(unlist(predictTBL[predictTBL$Management != "Fished", paste0("RealizedGainsSample_", i)])))
}

#Create a long dataframe with all of the samples
manageRealGainsSamplesDF <- rbind(data.frame(RealizedGains = restRG,
                                             Management = "Restricted"),
                                  data.frame(RealizedGains = unfishedLowRG,
                                             Management = "UnfishedLow"),
                                  data.frame(RealizedGains = unfishedHighSmallNewRG,
                                             Management = "UnfishedHighSmallNew"),
                                  data.frame(RealizedGains = unfishedHighBigOldRG,
                                             Management = "UnfishedHighBigOld")) %>% 
  mutate(RealizedGainsLogPos = ifelse(test = RealizedGains <= zeroLogTrans, yes = zeroLogTrans, no = RealizedGains),
         Management = factor(Management, levels = c("UnfishedHighBigOld", "UnfishedHighSmallNew", "UnfishedLow", "Restricted")))

###Get the median and quantiles for each protection category, for all high compliance MPAs, and all sites
manageRealGainQuant <- tibble(Management = c("Restricted", "UnfishedLow", "UnfishedHighSmallNew", "UnfishedHighBigOld", "AllHigh", "All")) %>% 
  mutate(Median = NA,
         Quant025 = NA,
         Quant25 = NA,
         Quant75 = NA,
         Quant975 = NA,
         ypos = NA)

manageRealGainQuant$Median[manageRealGainQuant$Management == "Restricted"] <- median(restRG)
manageRealGainQuant$Median[manageRealGainQuant$Management == "UnfishedLow"] <- median(unfishedLowRG)
manageRealGainQuant$Median[manageRealGainQuant$Management == "UnfishedHighSmallNew"] <- median(unfishedHighSmallNewRG)
manageRealGainQuant$Median[manageRealGainQuant$Management == "UnfishedHighBigOld"] <- median(unfishedHighBigOldRG)
manageRealGainQuant$Median[manageRealGainQuant$Management == "AllHigh"] <- median(c(unfishedHighBigOldRG, unfishedHighSmallNewRG))
manageRealGainQuant$Median[manageRealGainQuant$Management == "All"] <- median(allRG)

manageRealGainQuant$Quant025[manageRealGainQuant$Management == "Restricted"] <- unname(quantile(x = restRG, probs = c(0.025)))
manageRealGainQuant$Quant025[manageRealGainQuant$Management == "UnfishedLow"] <- unname(quantile(x = unfishedLowRG, probs = c(0.025)))
manageRealGainQuant$Quant025[manageRealGainQuant$Management == "UnfishedHighSmallNew"] <- unname(quantile(x = unfishedHighSmallNewRG, probs = c(0.025)))
manageRealGainQuant$Quant025[manageRealGainQuant$Management == "UnfishedHighBigOld"] <- unname(quantile(x = unfishedHighBigOldRG, probs = c(0.025)))
manageRealGainQuant$Quant025[manageRealGainQuant$Management == "AllHigh"] <- unname(quantile(x = c(unfishedHighBigOldRG, unfishedHighSmallNewRG), probs = c(0.025)))
manageRealGainQuant$Quant025[manageRealGainQuant$Management == "All"] <- unname(quantile(x = allRG, probs = c(0.025)))

manageRealGainQuant$Quant25[manageRealGainQuant$Management == "Restricted"] <- unname(quantile(x = restRG, probs = c(0.25)))
manageRealGainQuant$Quant25[manageRealGainQuant$Management == "UnfishedLow"] <- unname(quantile(x = unfishedLowRG, probs = c(0.25)))
manageRealGainQuant$Quant25[manageRealGainQuant$Management == "UnfishedHighSmallNew"] <- unname(quantile(x = unfishedHighSmallNewRG, probs = c(0.25)))
manageRealGainQuant$Quant25[manageRealGainQuant$Management == "UnfishedHighBigOld"] <- unname(quantile(x = unfishedHighBigOldRG, probs = c(0.25)))
manageRealGainQuant$Quant25[manageRealGainQuant$Management == "AllHigh"] <- unname(quantile(x = c(unfishedHighBigOldRG, unfishedHighSmallNewRG), probs = c(0.25)))
manageRealGainQuant$Quant25[manageRealGainQuant$Management == "All"] <- unname(quantile(x = allRG, probs = c(0.25)))

manageRealGainQuant$Quant75[manageRealGainQuant$Management == "Restricted"] <- unname(quantile(x = restRG, probs = c(0.75)))
manageRealGainQuant$Quant75[manageRealGainQuant$Management == "UnfishedLow"] <- unname(quantile(x = unfishedLowRG, probs = c(0.75)))
manageRealGainQuant$Quant75[manageRealGainQuant$Management == "UnfishedHighSmallNew"] <- unname(quantile(x = unfishedHighSmallNewRG, probs = c(0.75)))
manageRealGainQuant$Quant75[manageRealGainQuant$Management == "UnfishedHighBigOld"] <- unname(quantile(x = unfishedHighBigOldRG, probs = c(0.75)))
manageRealGainQuant$Quant75[manageRealGainQuant$Management == "AllHigh"] <- unname(quantile(x = c(unfishedHighBigOldRG, unfishedHighSmallNewRG), probs = c(0.75)))
manageRealGainQuant$Quant75[manageRealGainQuant$Management == "All"] <- unname(quantile(x = allRG, probs = c(0.75)))

manageRealGainQuant$Quant975[manageRealGainQuant$Management == "Restricted"] <- unname(quantile(x = restRG, probs = c(0.975)))
manageRealGainQuant$Quant975[manageRealGainQuant$Management == "UnfishedLow"] <- unname(quantile(x = unfishedLowRG, probs = c(0.975)))
manageRealGainQuant$Quant975[manageRealGainQuant$Management == "UnfishedHighSmallNew"] <- unname(quantile(x = unfishedHighSmallNewRG, probs = c(0.975)))
manageRealGainQuant$Quant975[manageRealGainQuant$Management == "UnfishedHighBigOld"] <- unname(quantile(x = unfishedHighBigOldRG, probs = c(0.975)))
manageRealGainQuant$Quant975[manageRealGainQuant$Management == "AllHigh"] <- unname(quantile(x = c(unfishedHighBigOldRG, unfishedHighSmallNewRG), probs = c(0.975)))
manageRealGainQuant$Quant975[manageRealGainQuant$Management == "All"] <- unname(quantile(x = allRG, probs = c(0.975)))

manageRealGainQuant$ypos[manageRealGainQuant$Management == "Restricted"] <- -100000
manageRealGainQuant$ypos[manageRealGainQuant$Management == "UnfishedLow"] <- -200000
manageRealGainQuant$ypos[manageRealGainQuant$Management == "UnfishedHighSmallNew"] <- -300000
manageRealGainQuant$ypos[manageRealGainQuant$Management == "UnfishedHighBigOld"] <- -400000

##Run a Kruskal-Wallis test with pairwise wilcox post tests
kwTestTBL <- manageRealGainsSamplesDF %>% 
  mutate(Management = gsub(pattern = "SmallNew|BigOld", replacement = "", x = as.character(Management)))

kruskal.test(RealizedGains ~ Management, data = kwTestTBL) #p << 0.001

pairwise.wilcox.test(kwTestTBL$RealizedGains, kwTestTBL$Management,
                     p.adjust.method = "BH") 

##Use a specific wilcox test to compare realized gains in all high compliance sites with the restricted sites
wilcox.test(RealizedGains ~ Management,
            data = manageRealGainsSamplesDF %>%
              filter(Management %in% c("Restricted", "UnfishedHighSmallNew", "UnfishedHighBigOld")) %>%
              mutate(Management = ifelse(Management == "Restricted", "Restricted", "UnfishedHigh")),
            conf.int = T, correct = F)

sumCounts <- numSamples*nrow(predictTBL)

freqBreaks <- c(0, 0.2*sumCounts, 0.4*sumCounts, 0.6*sumCounts)
freqLabels <- c(0, 0.2, 0.4, 0.6)

###>> Fig 2a plot ####
Fig2a_RealizedGainsStackedDistByManagePlot <- ggplot() + 
  geom_density(data = manageRealGainsSamplesDF,
               position = "stack",
               aes(x = RealizedGainsLogPos, y = after_stat(count), fill = Management)) +
  geom_errorbarh(data = manageRealGainQuant %>% filter(!Management %in% c("AllHigh", "All")) %>% mutate(Quant025 = ifelse(test = Quant025 < zeroLogTrans, yes = zeroLogTrans, no = Quant025),
                                                                                                        Quant975 = ifelse(test = Quant975 < zeroLogTrans, yes = zeroLogTrans, no = Quant975)),
                 aes(xmin = Quant025, xmax = Quant975, y = ypos, color = Management), linewidth = 1, height = 0) +
  geom_errorbarh(data = manageRealGainQuant %>% filter(!Management %in% c("AllHigh", "All")) %>% mutate(Quant25 = ifelse(test = Quant25 < zeroLogTrans, yes = zeroLogTrans, no = Quant25),
                                                                                                        Quant75 = ifelse(test = Quant75 < zeroLogTrans, yes = zeroLogTrans, no = Quant75)),
                 aes(xmin = Quant25, xmax = Quant75, y = ypos, color = Management), linewidth = 2, height = 0) +
  geom_point(data = manageRealGainQuant %>% filter(!Management %in% c("AllHigh", "All")) %>% mutate(Median = ifelse(test = Median < zeroLogTrans, yes = zeroLogTrans, no = Median)),
             aes(x = Median, y = ypos, fill = Management), shape = 21, size = 3) +
  geom_hline(yintercept = 0) +
  scale_x_continuous(trans = "log10", 
                     expand = c(0, 0), 
                     breaks = gainAmount_breaks,
                     labels = gainAmount_labels,
                     name = "Realized fish biomass gains (kg/ha)") +
  scale_y_continuous(name = "Relative frequency", labels = freqLabels, breaks = freqBreaks) +
  guides(fill = "none", color = "none") +
  scale_fill_manual(values = manageColors,
                    labels = manageLabels) +
  scale_color_manual(values = manageColors,
                     labels = manageLabels) +
  theme_bw() +
  annotation_logticks(sides = "b", outside = T) +
  coord_cartesian(clip = "off") +
  theme(axis.line = element_line(color = 'black', size = 1),
        axis.text.x = element_text(margin = margin(t = 8)),
        legend.position = c(0.6, 0.8),
        legend.title = element_blank(),
        plot.background = element_blank() ,
        panel.grid.major = element_blank() ,
        panel.grid.minor = element_blank() ,
        panel.border = element_blank() ,
        panel.background = element_blank(),
        #axis.text.y = element_blank(),
        #axis.ticks.y = element_blank(),
        aspect.ratio = 1) 

#Run again to get the legend
Fig2a_LegendPlot <- ggplot() + 
  geom_density(data = manageRealGainsSamplesDF,
               position = "stack",
               aes(x = RealizedGainsLogPos, y = after_stat(count), fill = Management)) +
  geom_errorbarh(data = manageRealGainQuant %>% filter(!Management %in% c("AllHigh", "All")) %>% mutate(Quant025 = ifelse(test = Quant025 < zeroLogTrans, yes = zeroLogTrans, no = Quant025),
                                                                                                        Quant975 = ifelse(test = Quant975 < zeroLogTrans, yes = zeroLogTrans, no = Quant975)),
                 aes(xmin = Quant025, xmax = Quant975, y = ypos, color = Management), linewidth = 1, height = 0) +
  geom_errorbarh(data = manageRealGainQuant %>% filter(!Management %in% c("AllHigh", "All")) %>% mutate(Quant25 = ifelse(test = Quant25 < zeroLogTrans, yes = zeroLogTrans, no = Quant25),
                                                                                                        Quant75 = ifelse(test = Quant75 < zeroLogTrans, yes = zeroLogTrans, no = Quant75)),
                 aes(xmin = Quant25, xmax = Quant75, y = ypos, color = Management), linewidth = 2, height = 0) +
  geom_point(data = manageRealGainQuant %>% filter(!Management %in% c("AllHigh", "All")) %>% mutate(Median = ifelse(test = Median < zeroLogTrans, yes = zeroLogTrans, no = Median)),
             aes(x = Median, y = ypos, fill = Management), shape = 21, size = 3) +
  geom_hline(yintercept = 0) +
  scale_x_continuous(trans = "log10", 
                     expand = c(0, 0), 
                     breaks = gainAmount_breaks,
                     labels = gainAmount_labels,
                     name = "Realized fish biomass gains (kg/ha)") +
  scale_y_continuous(name = "Frequency") +
  guides(fill = "none", color = "none") +
  scale_fill_manual(name = NULL,
                    values = manageColors,
                    labels = manageLabels) +
  scale_color_manual(name = NULL,
                     values = manageColors,
                     labels = manageLabels) +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 6), order = 1, reverse = T),
         size = "none") +
  theme_bw() +
  annotation_logticks(sides = "b", outside = T) +
  coord_cartesian(clip = "off") +
  theme(axis.line = element_line(color = 'black', size = 1),
        axis.text.x = element_text(margin = margin(t = 8)),
        legend.box.spacing = unit(1, 'cm'),
        legend.position = "top",
        legend.box = "vertical",
        legend.key = element_rect(fill = NA),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-20, -20, -20, -20),
        legend.text.align = 0.5,
        plot.background = element_blank() ,
        panel.grid.major = element_blank() ,
        panel.grid.minor = element_blank() ,
        panel.border = element_blank() ,
        panel.background = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        aspect.ratio = 1) 

Fig2a_Legend <- get_legend(Fig2a_LegendPlot)

######      > 2b. Realized coral reef fish biomass gains along a gradient of nearest market gravity ####
#Get predictions from each of the best models across the full gradient of gravity
gravGradientUnfishedBigOldTBL <- tibble(logGrav_NearMarket = seq(from = min(predictTBL$logGrav_NearMarket),
                                                                 to = max(predictTBL$logGrav_NearMarket),
                                                                 length.out = 1000)) %>% 
  mutate(Habitat = "Slope",
         Depth = "4-10m",
         CensusMethod = "Belt transect",
         Management = "UnfishedHighBigOld",
         logDHW_max_2yr = mean(predictTBL$logDHW_max_2yr),
         logChlA_mean_2yr = mean(predictTBL$logChlA_mean_2yr),
         logChlA_max_2yr = mean(predictTBL$logChlA_max_2yr),
         PAR_mean_2yr = mean(predictTBL$PAR_mean_2yr),
         PAR_sd_2yr = mean(predictTBL$PAR_sd_2yr),
         PAR_skewness_2yr = mean(predictTBL$PAR_skewness_2yr),
         logwave_mean = mean(predictTBL$logwave_mean),
         logGrav_NearPop = mean(predictTBL$logGrav_NearPop),
         SST_mean_2yr = mean(predictTBL$SST_mean_2yr),
         SST_max_2yr = mean(predictTBL$SST_max_2yr),
         SST_skewness_2yr = mean(predictTBL$SST_skewness_2yr),
         logPAR_kurtosis_2yr = mean(predictTBL$logPAR_kurtosis_2yr))

gravGradientfishedLessThan50kmTBL <- gravGradientUnfishedBigOldTBL %>% 
  mutate(Management = "Fished")

gravGradientPredsTBL <- tibble(logGrav_NearMarket = as.numeric(NA),
                               ModelNum = as.numeric(NA),
                               UnfishedBigOldPredMean = as.numeric(NA),
                               FishedPredMean = as.numeric(NA)) %>% 
  filter(!is.na(ModelNum))

i = 1
for(i in 1:nrow(bestSpammResTBL)) {
  message("Started model ", i, " of ", nrow(bestSpammResTBL))
  iterSpaMM <- readRDS(file = paste0(modelDir, bestSpammResTBL$ModelFilename[i]))
  
  iterGravGradientPredsTBL <- gravGradientUnfishedBigOldTBL %>% 
    select(logGrav_NearMarket) %>% 
    distinct() %>% 
    mutate(ModelNum = i,
           UnfishedBigOldLogPredMean = unname(predict(object = iterSpaMM,
                                                   type = "link",
                                                   binding = NA,
                                                   newdata = gravGradientUnfishedBigOldTBL,
                                                   re.form = NA)),
           FishedLogPredMean = unname(predict(object = iterSpaMM,
                                           type = "link",
                                           binding = NA,
                                           newdata = gravGradientfishedLessThan50kmTBL,
                                           re.form = NA)))
  gravGradientPredsTBL <- rbind(gravGradientPredsTBL, iterGravGradientPredsTBL)
}

gravGradientPredsTBL <- gravGradientPredsTBL %>% 
  mutate(RealizedGainsMean = exp(UnfishedBigOldLogPredMean) - exp(FishedLogPredMean),
         ModelNum = as.factor(ModelNum))

#Get the average biomasses for each gravity
gravGradientAllModelsTBL <- gravGradientPredsTBL %>% 
  group_by(logGrav_NearMarket) %>% 
  summarise(UnfishedBigOldLogPredMean = mean(UnfishedBigOldLogPredMean),
            FishedLogPredMean = mean(FishedLogPredMean)) %>% 
  mutate(RealizedGainsMean = exp(UnfishedBigOldLogPredMean) - exp(FishedLogPredMean))

###>> Fig 2b plot ####
Fig2b_RealizedGainsVsNearMarketGravMargSepModelsPlot <- ggplot(data = gravGradientPredsTBL,
                                                                      aes(x = logGrav_NearMarket, y = RealizedGainsMean)) +
  geom_line(color = manageColors["UnfishedHighBigOld"], linewidth = 1, linetype = "dotted", aes(group = ModelNum)) +
  geom_line(data = gravGradientAllModelsTBL, aes(x = logGrav_NearMarket, y = RealizedGainsMean), linewidth = 2, color = manageColors["UnfishedHighBigOld"]) +
  scale_x_continuous(name = "Nearest market gravity",
                     breaks = c(log(0.0001), log(0.01), log(1), log(100), log(10000)),
                     labels = c(0.0001, 0.01, 1, 100, 10000)) +
  scale_y_continuous(name = "Realized fish biomass gains (kg/ha)") +
  guides(linetype = guide_legend(title = "Model")) +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.subtitle = element_text(hjust=0.5),
        axis.line = element_line(color = 'black'),
        aspect.ratio = 1)

######      > Assemble figure 2 ####
Fig2_RealGainsStackedHistGradientGravity <- ggpubr::ggarrange(
  #First row with the cumulative plot and stacked histogram for the full protection gains
  ggpubr::ggarrange(Fig2a_RealizedGainsStackedDistByManagePlot,
                    Fig2b_RealizedGainsVsNearMarketGravMargSepModelsPlot,
                    labels = c("a", "b"), ncol = 2),
  #Second row with fill legend
  ggpubr::ggarrange(Fig2a_Legend),
  nrow = 2, heights = c(2.1,0.2)
)

#Plot as a tiff
ggsave(plot = Fig2_RealGainsStackedHistGradientGravity,
       filename = paste0(plotDir, "Caldwelletal_RealizedPotentialGains_Fig2_RealGainsStackedHistGradientGravity.tiff"),
       dpi = 1000, device = "tiff", bg = "white", width = 10, height = 5)

#Plot as a PDF
ggsave(plot = Fig2_RealGainsStackedHistGradientGravity,
       filename = paste0(plotDir, "Caldwelletal_RealizedPotentialGains_Fig2_RealGainsStackedHistGradientGravity.pdf"),
       dpi = 1000, device = "pdf", bg = "white", width = 10, height = 5)

######    Fig 3. Potential gains ####
#Get the median and quantiles for all the unfished scenario to add to paper
allUnfishedPredVect <- as.vector(t(select(predictTBL, starts_with("UnfishedPredSample_"))))
median(allUnfishedPredVect) #956
quantile(allUnfishedPredVect, probs = c(0.025, 0.975)) #404 to 2766 kg/ha

######      > 3a. Cumulative potential and realized gains (% of status quo - potential gains positive, realized gains negative) ####
cumGainsTBL <- bind_rows(predictTBL %>%
                           dplyr::filter(Management %in% c("UnfishedHighBigOld", "UnfishedHighSmallNew")) %>% 
                           arrange(desc(Management), desc(RealizedGainMed)),
                         predictTBL %>% 
                           dplyr::filter(!Management %in% c("UnfishedHighBigOld", "UnfishedHighSmallNew")) %>% 
                           arrange(desc(PotentialGainMed))) %>% 
  mutate(PercentageSites = c(1:nrow(predictTBL))/nrow(predictTBL)*100) %>% 
  select(Management, PercentageSites, contains("PotentialGainsPercStatusQuoSample_"), contains("RealizedGainsPercStatusQuoSample_")) 

cumPotGainsPercStatusQuoQuantiles <- select(cumGainsTBL, starts_with("PotentialGainsPercStatusQuoSample_")) %>%
  mutate(across(everything(), ~cumsum(.))) %>%
  rename_with(~paste0(., "_cumulative"), everything()) %>%
  rowwise() %>%
  mutate(
    CumPotGainsPercTotalStatusQuoMed = median(c_across(starts_with("PotentialGainsPercStatusQuoSample_")), na.rm = TRUE),
    CumPotGainsPercTotalStatusQuoQuant_025 = quantile(c_across(starts_with("PotentialGainsPercStatusQuoSample_")), probs = 0.025, na.rm = TRUE),
    CumPotGainsPercTotalStatusQuoQuant_975 = quantile(c_across(starts_with("PotentialGainsPercStatusQuoSample_")), probs = 0.975, na.rm = TRUE)) %>% 
  select(CumPotGainsPercTotalStatusQuoMed, CumPotGainsPercTotalStatusQuoQuant_025, CumPotGainsPercTotalStatusQuoQuant_975)

cumRealGainsPercStatusQuoQuantiles <- select(cumGainsTBL, starts_with("RealizedGainsPercStatusQuoSample_")) %>%
  mutate(across(everything(), ~cumsum(.))) %>%
  rename_with(~paste0(., "_cumulative"), everything()) %>%
  rowwise() %>%
  mutate(
    CumRealGainsPercTotalStatusQuoMed = median(c_across(starts_with("RealizedGainsPercStatusQuoSample_")), na.rm = TRUE)*-1,
    CumRealGainsPercTotalStatusQuoQuant_025 = quantile(c_across(starts_with("RealizedGainsPercStatusQuoSample_")), probs = 0.025, na.rm = TRUE)*-1,
    CumRealGainsPercTotalStatusQuoQuant_975 = quantile(c_across(starts_with("RealizedGainsPercStatusQuoSample_")), probs = 0.975, na.rm = TRUE)*-1) %>% 
  select(CumRealGainsPercTotalStatusQuoMed, CumRealGainsPercTotalStatusQuoQuant_025, CumRealGainsPercTotalStatusQuoQuant_975)

cumGainsTBL <- bind_cols(cumGainsTBL, cumPotGainsPercStatusQuoQuantiles, cumRealGainsPercStatusQuoQuantiles)

#Calculate cumulative potential gains for randomly selected sites - to compare with ordered
highCompGainsTBL <- cumGainsTBL %>% 
  dplyr::filter(Management %in% c("UnfishedHighSmallNew", "UnfishedHighBigOld"))

noHighCompGainsTBL <- cumGainsTBL %>% 
  dplyr::filter(!Management %in% c("UnfishedHighSmallNew", "UnfishedHighBigOld"))

randCumPotGainsFilename <- paste0(dataDir, "Caldwelletal_RealizedPotentialGains_RandomCumulativePotentialGains.rds")
if(file.exists(randCumPotGainsFilename)) {
  randCumPotentialGainsTBL <- readRDS(file = randCumPotGainsFilename)
} else {
  randCumPotentialGainsTBL <- data.frame(matrix(nrow = nrow(cumGainsTBL), ncol = numSamples)) 
  
  iterNum <- 0
  i <- 1
  for (i in 1:numSamples) {
    iterNum <- iterNum + 1
    
    randNoHighRows <- sample(nrow(noHighCompGainsTBL))
    randNoHighTempTBL <- noHighCompGainsTBL[randNoHighRows,]
    randAllTempTBL <- bind_rows(highCompGainsTBL, randNoHighTempTBL)
    
    potentialGainsColname <- paste0("PotentialGainsPercStatusQuoSample_",i)
    
    randCumPotentialGainsTBL[,iterNum] <- cumsum(randAllTempTBL[,potentialGainsColname])
    message("Random cumulative gains = ", max(randCumPotentialGainsTBL[,iterNum]))
    
  }

  #Save this dataset
  saveRDS(object = randCumPotentialGainsTBL, file = randCumPotGainsFilename)
}

randCumPotentialGainsQuantilesTBL <- randCumPotentialGainsTBL %>%
  rowwise() %>%
  mutate(
    RandCumPotGainsCompStatusQuoMed = median(c_across(starts_with("X")), na.rm = TRUE),
    RandCumPotGainsCompStatusQuoQuant_025 = quantile(c_across(starts_with("X")), probs = 0.025, na.rm = TRUE),
    RandCumPotGainsCompStatusQuoQuant_975 = quantile(c_across(starts_with("X")), probs = 0.975, na.rm = TRUE)
  ) %>% 
  select(RandCumPotGainsCompStatusQuoMed, RandCumPotGainsCompStatusQuoQuant_025, RandCumPotGainsCompStatusQuoQuant_975) %>% 
  bind_cols(cumGainsTBL %>% select(PercentageSites))

##Find out the median and 95% quantiles for expected potential gains at ~30% high compliance
#Get the percentage closest to 30%
closestPerc30 <- randCumPotentialGainsQuantilesTBL$PercentageSites[which.min(abs(randCumPotentialGainsQuantilesTBL$PercentageSites - 30))]

#Maximizing gains
maxPotGains30PercMed <- cumGainsTBL$CumPotGainsPercTotalStatusQuoMed[cumGainsTBL$PercentageSites == closestPerc30] #28.2%
maxPotGains30PercQuant025 <- cumGainsTBL$CumPotGainsPercTotalStatusQuoQuant_025[cumGainsTBL$PercentageSites == closestPerc30] #26.1
maxPotGains30PercQuant975 <- cumGainsTBL$CumPotGainsPercTotalStatusQuoQuant_975[cumGainsTBL$PercentageSites == closestPerc30] #30.4

#Random
randPotGains30PercMed <- randCumPotentialGainsQuantilesTBL$RandCumPotGainsCompStatusQuoMed[randCumPotentialGainsQuantilesTBL$PercentageSites == closestPerc30] #16.6
randPotGains30PercQuant025 <- randCumPotentialGainsQuantilesTBL$RandCumPotGainsCompStatusQuoQuant_025[randCumPotentialGainsQuantilesTBL$PercentageSites == closestPerc30] #15.2
randPotGains30PercQuant975 <- randCumPotentialGainsQuantilesTBL$RandCumPotGainsCompStatusQuoQuant_975[randCumPotentialGainsQuantilesTBL$PercentageSites == closestPerc30] #18.1

####>> Fig 3a plot ####
Fig3a_CumPotRealGainsVsPerSites <- ggplot(data = cumGainsTBL,
                                          aes(x = PercentageSites)) +
  #Add a ribbon for the range of values for random siting
  geom_ribbon(data = randCumPotentialGainsQuantilesTBL %>%
                filter(PercentageSites > max(highCompGainsTBL$PercentageSites)),
              aes(ymin = RandCumPotGainsCompStatusQuoQuant_025,
                  ymax = RandCumPotGainsCompStatusQuoQuant_975),
              color = "grey", fill = "grey",
              show.legend = F) +
  #Add a dashed line for the mean random siting
  geom_line(data = randCumPotentialGainsQuantilesTBL %>%
              filter(PercentageSites > max(highCompGainsTBL$PercentageSites)),
            aes(y = RandCumPotGainsCompStatusQuoMed),
            linetype = "dashed", size = 1, color = "black") +
  #Rectangle showing realized and potential gains
  geom_rect(aes(xmin = 100, xmax = 105, ymin = min(CumRealGainsPercTotalStatusQuoQuant_975), ymax = 0),
            color = "black", fill = gainColors["RealizedGains"]) +
  geom_rect(aes(xmin = 100, xmax = 105, ymin = 0, ymax = max(CumPotGainsPercTotalStatusQuoQuant_975)),
            color = "black", fill = gainColors["PotentialGains"]) +
  geom_segment(x = 100, xend = 105,
               aes(y = max(CumPotGainsPercTotalStatusQuoMed),
                   yend = max(CumPotGainsPercTotalStatusQuoMed)),
               color = manageColors["UnfishedHighBigOld"], size = 1) +
  geom_segment(x = 100, xend = 105,
               aes(y = min(CumRealGainsPercTotalStatusQuoMed),
                   yend = min(CumRealGainsPercTotalStatusQuoMed)),
               color = manageColors["Fished"], size = 1) +
  geom_segment(x = 100, xend = 105, y = 0, yend = 0,
               color = "black", size = 2) +
  geom_ribbon(aes(ymin = CumRealGainsPercTotalStatusQuoQuant_025,
                  ymax = CumRealGainsPercTotalStatusQuoQuant_975),
              fill = manageColors["Fished"],
              alpha = 0.5) +
  geom_line(aes(y = CumRealGainsPercTotalStatusQuoMed),
            color = manageColors["Fished"],
            size = 1, linetype = "dashed") +
  geom_ribbon(aes(ymin = CumPotGainsPercTotalStatusQuoQuant_025,
                  ymax = CumPotGainsPercTotalStatusQuoQuant_975),
              fill = manageColors["UnfishedHighBigOld"],
              alpha = 0.5) +
  geom_line(aes(y = CumPotGainsPercTotalStatusQuoMed),
            color = manageColors["UnfishedHighBigOld"],
            size = 1, linetype = "dashed") +
  geom_point(aes(color = Management), y = 0, size = 2) +
  guides(color = "none", fill = "none") +
  scale_x_continuous(name = "% of sites",
                     limits = c(0, 105),
                     expand = c(0, 0),
                     breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) +
  scale_y_continuous(name = "Cumulative fish biomass gains\n(% of status quo biomass)",
                     expand = c(0,0),
                     breaks = c(-30, -20, -10, 0, 10, 20, 30, 40, 50, 60, 70, 80),
                     sec.axis = dup_axis(breaks = c(max(cumGainsTBL$CumPotGainsPercTotalStatusQuoMed),
                                                    0,
                                                    min(cumGainsTBL$CumRealGainsPercTotalStatusQuoMed)),
                                         labels = c(paste0("Full protection (+", round(max(cumGainsTBL$CumPotGainsPercTotalStatusQuoMed)), "%)"),
                                                    "Status quo",
                                                    paste0("All fished (", round(min(cumGainsTBL$CumRealGainsPercTotalStatusQuoMed)),"%)")),
                                         name = NULL)) +
  scale_color_manual(values = manageColors,
                     labels = manageLabels) +
  scale_fill_manual(values = manageColors,
                    labels = manageLabels) +
  #Add dotted lines for maximum and random potential gains at 30% protection
  geom_segment(x = closestPerc30, xend = closestPerc30,
               y = 0, yend = maxPotGains30PercMed,
               color = "black", linetype = "dotted") +
  geom_segment(x = 0, xend = closestPerc30,
               y = maxPotGains30PercMed, yend = maxPotGains30PercMed,
               color = "black", linetype = "dotted") +
  geom_segment(x = 0, xend = closestPerc30,
               y = randPotGains30PercMed, yend = randPotGains30PercMed,
               color = "black", linetype = "dotted") +
  theme_bw() +
  theme(plot.background = element_blank() ,
        panel.grid.major = element_blank() ,
        panel.grid.minor = element_blank() ,
        panel.background = element_blank(),
        axis.line = element_line(color = 'black'),
        aspect.ratio = 1) +
  #Add text for realized and potential gains
  geom_text(x = 102.5,
            y = max(cumGainsTBL$CumPotGainsPercTotalStatusQuoMed)/2,
            label = "Potential",
            color = "white",
            hjust = 0.5,
            vjust = 0.5,
            angle = 270) +
  geom_text(x = 102.5,
            y = min(cumGainsTBL$CumRealGainsPercTotalStatusQuoMed)/2,
            label = "Realized",
            color = "black",
            hjust = 0.5,
            vjust = 0.5,
            angle = 270)

######      > 3b. Map of potential gains ####
###Prepare data
potGainMapTBL <- predictTBL %>%
  select(Long_DD, Lat_DD, Management, PotentialGainMed) %>% 
  mutate(Long_DD = ifelse(test = Long_DD < -25, yes = Long_DD + 360, no = Long_DD),
         PotentialGainsLogPos = ifelse(test = PotentialGainMed < zeroLogTrans, yes = zeroLogTrans, no = PotentialGainMed),
         Management = factor(Management, levels = c("UnfishedHighBigOld", "UnfishedHighSmallNew", "UnfishedLow", "Restricted", "Fished"))) %>% 
  arrange(desc(PotentialGainsLogPos), desc(Management))

####>> Fig 3b map #### 
Fig3b_PotGainsMap <- ggplot() + 
  geom_polygon(data = mapWorld, aes(x = long, y = lat, group = group), fill = "grey50", color = "grey50") +
  coord_fixed(ratio = 1.3, xlim = c(40, 320), ylim = c(-25, 25)) +
  geom_jitter(data = potGainMapTBL,
              aes(x = Long_DD, y = Lat_DD, fill = Management, size = PotentialGainsLogPos),
              colour = "black", pch = 21, width = 1, height = 1) +
  scale_fill_manual(name = "Current (status quo) management:",
                    values = manageColors,
                    labels = manageLabels) +
  scale_size_continuous(breaks = mapGain_breaks, 
                        labels = mapGain_labels, 
                        range = c(1,8), 
                        limits = gainAmount_limits) +
  guides(fill = "none",
         size = guide_legend(title = "Potential fish biomass gains (kg/ha):", order = 2)) +
  geom_hline(yintercept = 23.43695, lty = 2) +
  geom_hline(yintercept = -23.43695, lty = 2) + 
  scale_x_continuous("", breaks = c(30, 80, 130, 180, 230, 280, 330), labels=c(-150, -100, -50, 0, 50, 100, 150)) +
  scale_y_continuous("", breaks = c(-20, -10, 0, 10, 20)) +
  theme(plot.subtitle = element_text(face = "italic"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        #legend.box.spacing = unit(1, 'cm'),
        legend.title.align = 0.5,
        legend.background = element_rect(fill = "transparent"),
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.key = element_blank(),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-20, -20, -20, -20),
        legend.text.align = 0.5)

######      > Assemble figure 3 ####
Fig3_PotGainsCumBiomassMap <- ggpubr::ggarrange(
  #First row with the cumulative plot and stacked histogram for the full protection gains
  ggpubr::ggarrange(Fig3a_CumPotRealGainsVsPerSites, labels = c("a")),
  #Second row with fill legend
  ProtectionCategory_FillLegend,
  #Third row with the map
  ggpubr::ggarrange(Fig3b_PotGainsMap, labels = c("b")),
  nrow = 3, heights = c(2.1,0.1,1.5)
)

#Plot as a tiff
ggsave(plot = Fig3_PotGainsCumBiomassMap,
       filename = paste0(plotDir, "Caldwelletal_RealizedPotentialGains_Fig3_PotGainsCumGainsMap.tiff"),
       dpi = 1000, device = "tiff", bg = "white", width = 11, height = 11)

#Plot as a PDF
ggsave(plot = Fig3_PotGainsCumBiomassMap,
       filename = paste0(plotDir, "Caldwelletal_RealizedPotentialGains_Fig3_PotGainsCumGainsMap.pdf"),
       dpi = 1000, device = "pdf", bg = "white", width = 11, height = 11)

####Calculate the amount of gains if all fished sites were in restricted fishing areas ####
restGainsFishedPercStatusQuoSums <- predictTBL %>% 
  filter(Management %in% c("Fished")) %>% 
  select(starts_with("RestGainsPercStatusQuoSample_"))

restGainsFishedPercStatusQuoSums[restGainsFishedPercStatusQuoSums < 0] <- 0 #truncate at zero to make comparable with potential gains

restGainsFishedPercStatusQuoSums <- restGainsFishedPercStatusQuoSums %>%
  summarize(across(everything(), ~sum(.)))

median(as.numeric(restGainsFishedPercStatusQuoSums)) #10.5% --> potential gains if all fished sites were restricted fishing
quantile(as.numeric(restGainsFishedPercStatusQuoSums), probs = 0.025) #9.9
quantile(as.numeric(restGainsFishedPercStatusQuoSums), probs = 0.975) #11.1

####SUPPLEMENTAL FIGURES ####
######    Fig S1. Model data and results: Status quo biomass (map and distribution plot) and effect sizes ####
#Get the range of median predictions for status quo to add to the paper
range(predictTBL$StatusQuoBiomassMed) #191 to 2387

######        S1a. Map the status quo biomasses for all sites (corrected by method, habitat, depth) ####
###Prepare data
statusQuoBiomassMapTBL <- predictTBL %>%
  select(Long_DD, Lat_DD, Management, StatusQuoBiomassMed) %>% 
  mutate(Long_DD = ifelse(test = Long_DD < -25, yes = Long_DD + 360, no = Long_DD),
         Management = factor(Management, levels = c("UnfishedHighBigOld", "UnfishedHighSmallNew", "UnfishedLow", "Restricted", "Fished"))) %>% 
  arrange(desc(StatusQuoBiomassMed))

#Create breaks
mapStatusQuoBiomass_breaks <- c(500, 1000, 1500, 2000)

###>> Fig S1a map ####
FigS1a_StatusQuoBiomassMapSize <- ggplot() +
  geom_polygon(data = mapWorld, aes(x = long, y = lat, group = group), fill = "grey50", color = "grey50") +
  coord_fixed(ratio = 1.3, xlim = c(40, 320), ylim = c(-25, 25)) +
  geom_jitter(data = statusQuoBiomassMapTBL,
              aes(x = Long_DD, y = Lat_DD, fill = Management, size = StatusQuoBiomassMed), 
              colour = "black", pch = 21, width = 1, height = 1) +
  scale_fill_manual(name = NULL,
                    values = manageColors,
                    labels = manageLabels) +
  scale_size_continuous(name = "Status quo fish biomass (kg/ha):",
                        breaks = mapStatusQuoBiomass_breaks,
                        labels = mapStatusQuoBiomass_breaks,
                        range = c(1, 8)) +
  guides(fill = guide_legend(override.aes = list(size = 5))) +
  geom_hline(yintercept = 23.43695, lty = 2)+
  geom_hline(yintercept = -23.43695, lty = 2)+
  scale_x_continuous("", breaks = c(30, 80, 130, 180, 230, 280, 330), labels=c(-150, -100, -50, 0, 50, 100, 150)) +
  scale_y_continuous("", breaks = c(-20, -10, 0, 10, 20)) +
  theme(plot.subtitle = element_text(face = "italic"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.box = "vertical",
        legend.key = element_rect(fill = NA),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-20, -20, -20, -20),
        legend.text.align = 0.5)

######        S1b. Stacked frequency distribution of status quo biomass ####
#Make a dataset with the medians and quantiles for each management category (from the samples)
#Create a long dataframe with all of the status quo biomas samples from each management type
manageStatusQuoSamplesTBL <- bind_rows(tibble(StatusQuoBiomass = as.vector(t(predictTBL %>% 
                                                                               filter(Management == "Restricted") %>% 
                                                                               select(starts_with("StatusQuoPredSample_")))),
                                              Management = "Restricted"),
                                       tibble(StatusQuoBiomass = as.vector(t(predictTBL %>% 
                                                                               filter(Management == "UnfishedLow") %>% 
                                                                               select(starts_with("StatusQuoPredSample_")))),
                                              Management = "UnfishedLow"),
                                       tibble(StatusQuoBiomass = as.vector(t(predictTBL %>% 
                                                                               filter(Management == "UnfishedHighSmallNew") %>% 
                                                                               select(starts_with("StatusQuoPredSample_")))),
                                              Management = "UnfishedHighSmallNew"),
                                       tibble(StatusQuoBiomass = as.vector(t(predictTBL %>% 
                                                                               filter(Management == "UnfishedHighBigOld") %>% 
                                                                               select(starts_with("StatusQuoPredSample_")))),
                                              Management = "UnfishedHighBigOld"),
                                       tibble(StatusQuoBiomass = as.vector(t(predictTBL %>% 
                                                                               filter(Management == "Fished") %>% 
                                                                               select(starts_with("StatusQuoPredSample_")))),
                                              Management = "Fished")) %>% 
  mutate(Management = factor(Management, levels = c("UnfishedHighBigOld", "UnfishedHighSmallNew", "UnfishedLow", "Restricted", "Fished")))

###Get the median and quantiles for each protection category, for all high compliance MPAs, and all sites
manageStatusQuoQuantTBL <- manageStatusQuoSamplesTBL %>% 
  group_by(Management) %>% 
  summarise(Median = median(StatusQuoBiomass),
            Quant025 = unname(quantile(StatusQuoBiomass, probs = c(0.025))),
            Quant25 = unname(quantile(StatusQuoBiomass, probs = c(0.25))),
            Quant75 = unname(quantile(StatusQuoBiomass, probs = c(0.75))),
            Quant975 = unname(quantile(StatusQuoBiomass, probs = c(0.975)))) %>% 
  ungroup() %>% 
  mutate(ypos = case_match(Management,
                           "UnfishedHighBigOld" ~ -400000,
                           "UnfishedHighSmallNew" ~ -800000,
                           "UnfishedLow" ~ -1200000,
                           "Restricted" ~ -1600000,
                           "Fished" ~ -2000000))
nrow(manageStatusQuoSamplesTBL)
sumCounts <- numSamples*nrow(predictTBL)

freqBreaks <- c(0, 0.2*sumCounts, 0.4*sumCounts, 0.6*sumCounts, 0.8*sumCounts, sumCounts)*2
freqLabels <- c(0, 0.2, 0.4, 0.6, 0.8, 1)

###>> Fig S1b plot ####
FigS1b_StatusQuoBiomassStackedDistByManagePlot <- ggplot() + 
  geom_density(data = manageStatusQuoSamplesTBL,
               position = "stack",
               aes(x = StatusQuoBiomass, y = after_stat(count), fill = Management)) +
  scale_y_continuous(name = "Relative frequency", labels = freqLabels, breaks = freqBreaks) +
  geom_errorbarh(data = manageStatusQuoQuantTBL,
                 aes(xmin = Quant025, xmax = Quant975, y = ypos, color = Management), linewidth = 1, height = 0) +
  geom_errorbarh(data = manageStatusQuoQuantTBL,
                 aes(xmin = Quant25, xmax = Quant75, y = ypos, color = Management), linewidth = 2, height = 0) +
  geom_point(data = manageStatusQuoQuantTBL,
             aes(x = Median, y = ypos, fill = Management), shape = 21, size = 3) +
  geom_hline(yintercept = 0) +
  scale_x_continuous(trans = "log10", 
                     expand = c(0, 0), 
                     # breaks = gainAmount_breaks,
                     # labels = gainAmount_labels,
                     name = "Status quo fish biomass (kg/ha)") +
  guides(fill = "none", color = "none") +
  scale_fill_manual(values = manageColors,
                    labels = manageLabels) +
  scale_color_manual(values = manageColors,
                     labels = manageLabels) +
  theme_bw() +
  annotation_logticks(sides = "b", outside = T) +
  coord_cartesian(clip = "off") +
  theme(axis.line = element_line(color = 'black', size = 1),
        axis.text.x = element_text(margin = margin(t = 8)),
        legend.title = element_blank(),
        plot.background = element_blank() ,
        panel.grid.major = element_blank() ,
        panel.grid.minor = element_blank() ,
        panel.border = element_blank() ,
        panel.background = element_blank(),
        #axis.text.y = element_blank(),
        #axis.ticks.y = element_blank(),
        aspect.ratio = 1) 

######        S1c. Plot the effect sizes ####
##Extract the fixed effects coefficients, SE, and degrees of freedom for each of the best models
standBetaTableDF <- tibble(Covariate = as.character(NA),
                           Estimate = as.numeric(NA),
                           CondSE = as.numeric(NA),
                           ModelIter = as.numeric(NA)) %>% 
  filter(!is.na(Covariate))

i = 1
for(i in 1:nrow(bestSpammResTBL)) {
  message("Started model ", i, " of ", nrow(bestSpammResTBL))
  #Open the best model
  iterBestStandSpammModel <- readRDS(file = paste0(modelDir, bestSpammResTBL$ModelFilename[i]))
  
  iterStandBetaTableDF <- spaMM::summary.HLfit(iterBestStandSpammModel)$beta_table %>% 
    as.data.frame() %>% 
    rownames_to_column() %>% 
    rename(Covariate = rowname,
           CondSE = `Cond. SE`) %>% 
    select(-"t-value") %>% 
    mutate(ModelIter = i)
  
  standBetaTableDF <- rbind(standBetaTableDF, iterStandBetaTableDF)
}

#Calculate the means and propagate error
meanStandBetaTBL <- standBetaTableDF %>% 
  group_by(Covariate) %>% 
  summarise(MeanEstimate = mean(Estimate, na.rm = T),
            WithinEstPropSE = sqrt(mean(CondSE^2, na.rm = T)),
            AmongEstSE = std.error(Estimate),
            NumModels = length(Estimate)) %>% 
  ungroup() %>% 
  mutate(AmongEstSE = replace_na(AmongEstSE, 0),
         TotalPropSE = sqrt((WithinEstPropSE^2) + (AmongEstSE^2)),
         ConfInt95 = 1.96 * TotalPropSE,
         Upper95 = MeanEstimate + ConfInt95,
         Lower95 = MeanEstimate - ConfInt95)

meanStandBetaTBL <- meanStandBetaTBL %>% 
  mutate(Color = case_when(Lower95 < 0 & Upper95 < 0 ~ "red",
                           Lower95 > 0 & Upper95 > 0 ~ "blue",
                           TRUE ~ "gray"),
         Group = case_when(grepl(pattern = "CensusMethod", x = Covariate) ~ 1,
                           grepl(pattern = "Habitat", x = Covariate) ~ 2,
                           grepl(pattern = "Depth", x = Covariate) ~ 3,
                           grepl(pattern = "Management", x = Covariate) & !grepl(pattern = "Grav", x = Covariate) ~ 4,
                           !grepl(pattern = "Management", x = Covariate) & grepl(pattern = "Grav", x = Covariate) ~ 5,
                           grepl(pattern = "Management", x = Covariate) & grepl(pattern = "Grav", x = Covariate) ~ 6,
                           grepl(pattern = "PAR", x = Covariate) ~ 7,
                           grepl(pattern = "ChlA", x = Covariate) ~ 8,
                           grepl(pattern = "SST", x = Covariate) ~ 9,
                           grepl(pattern = "DHW", x = Covariate) ~ 9,
                           grepl(pattern = "wave", x = Covariate) ~ 10),
         Label = case_when(Covariate == "HabitatLagoon/Back reef" ~ "Reef lagoon",
                           Covariate == "HabitatCrest" ~ "Reef crest",
                           Covariate == "HabitatFlat" ~ "Reef flat",
                           Covariate == "Depth>10m" ~ "Depth (>10m)",
                           Covariate == "Depth0-4m" ~ "Depth (0-4m)",
                           Covariate == "CensusMethodDistance sampling" ~ "Distance sampling method",
                           Covariate == "CensusMethodPoint intercept" ~ "Point intercept method",
                           Covariate == "ManagementRestricted" ~ "Fishing restricted",
                           Covariate == "ManagementUnfishedLow" ~ "Low compliance MPA",
                           Covariate == "ManagementUnfishedHighBigOld" ~ "Big and old high compliance MPA",
                           Covariate == "ManagementUnfishedHighSmallNew" ~ "Small or new high compliance MPA",
                           Covariate == "scale(logGrav_NearMarket)" ~ "Nearest market gravity",
                           Covariate == "scale(logGrav_NearPop)" ~ "Nearest population gravity",
                           Covariate == "scale(logwave_mean)" ~ "wave energy mean",
                           Covariate == "scale(PAR_skewness_2yr)" ~ "PAR skewness",
                           Covariate == "scale(PAR_mean_2yr)" ~ "PAR mean",
                           Covariate == "scale(PAR_sd_2yr)" ~ "PAR sd",
                           Covariate == "scale(logPAR_kurtosis_2yr)" ~ "PAR kurtosis",
                           Covariate == "scale(SST_max_2yr)" ~ "SST maximum",
                           Covariate == "scale(SST_mean_2yr)" ~ "SST mean",
                           Covariate == "scale(SST_skewness_2yr)" ~ "SST skewness",
                           Covariate == "scale(logDHW_max_2yr)" ~ "DHW max",
                           Covariate == "scale(logChlA_max_2yr)" ~ "Chl-a max",
                           Covariate == "scale(logChlA_mean_2yr)" ~ "Chl-a mean",
                           Covariate == "scale(logGrav_NearMarket):ManagementRestricted" ~ "Fishing restricted\n x Nearest market gravity",
                           Covariate == "scale(logGrav_NearMarket):ManagementUnfishedLow" ~ "Low compliance MPA\n x Nearest market gravity",
                           Covariate == "scale(logGrav_NearMarket):ManagementUnfishedHighBigOld" ~ "Big and old high compliance MPA\n x Nearest market gravity",
                           Covariate == "scale(logGrav_NearMarket):ManagementUnfishedHighSmallNew" ~ "Small or new high compliance MPA\n x Nearest market gravity"),
         Label = ifelse(NumModels < 3, yes = paste0(Label, " (", NumModels, ")"), no = Label)) %>% 
  arrange(Group, Color, desc(MeanEstimate)) %>% 
  mutate(Covariate = factor(x = Covariate, levels = rev(Covariate)),
         Label = factor(x = Label, levels = rev(Label))) %>% 
  filter(!is.na(Group))

####Plot the effect sizes 
#Common limits
effSizeLims <- c(min(meanStandBetaTBL$Lower95), max(meanStandBetaTBL$Upper95))

manageGravStdEffPlot <- ggplot() +
  geom_pointrange(meanStandBetaTBL %>% filter(Group %in% c(4,5,6)), 
                  mapping = aes(x = Label,
                                y = MeanEstimate,
                                colour = Color,
                                ymin = Lower95,
                                ymax = Upper95)) +
  scale_colour_manual(values = c("blue" = "blue",
                                 "gray" = "gray",
                                 "red" = "red")) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  geom_hline(yintercept = 0, color = "dark gray") +
  geom_vline(xintercept = 4.5, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 6.5, linetype = "dashed", color = "gray") +
  scale_y_continuous(limits = effSizeLims) +
  xlab("") +
  ylab("") +
  coord_flip()

methodHabitatDepthEffPlot <- ggplot() +
  geom_pointrange(meanStandBetaTBL %>% filter(Group %in% c(1,2,3)), 
                  mapping = aes(x = Label,
                                y = MeanEstimate,
                                colour = Color,
                                ymin = Lower95,
                                ymax = Upper95)) +
  scale_colour_manual(values = c("blue" = "blue",
                                 "gray" = "gray",
                                 "red" = "red")) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  geom_hline(yintercept = 0, color = "dark gray") +
  geom_vline(xintercept = 2.5, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 5.5, linetype = "dashed", color = "gray") +
  scale_y_continuous(limits = effSizeLims) +
  xlab("") +
  ylab("") +
  coord_flip()

envEffPlot <- ggplot() +
  geom_pointrange(meanStandBetaTBL %>% filter(Group %in% c(7,8,9,10)), 
                  mapping = aes(x = Label,
                                y = MeanEstimate,
                                colour = Color,
                                ymin = Lower95,
                                ymax = Upper95)) +
  scale_colour_manual(values = c("blue" = "blue",
                                 "gray" = "gray",
                                 "red" = "red")) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  geom_hline(yintercept = 0, color = "dark gray") +
  geom_vline(xintercept = 1.5, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 5.5, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 7.5, linetype = "dashed", color = "gray") +
  scale_y_continuous(limits = effSizeLims) +
  xlab("") +
  ylab("") +
  coord_flip()


######        S1d. Combine into one figure (several versions) ####
###Combine the two effect size plots first then combine with the rest
###>> Fig S1c plot ####
FigS1c_allCoeffPlot <- egg::ggarrange(manageGravStdEffPlot,
                               methodHabitatDepthEffPlot,
                               envEffPlot,
                               ncol = 3,
                               bottom = text_grob(label = "Standardized effect size", hjust = -0.15),
                               labels = c("c", "", ""), label.args =  list(gp = grid::gpar(font = 2, cex = 1.2)))

####Above: Map with color for management and size for biomass (alone); Below: Stacked histogram beside two coefficient plots 
FigS1_MapSizeColorHistGradCoeff <- ggpubr::ggarrange(
  ggpubr::ggarrange(FigS1a_StatusQuoBiomassMapSize, labels = c("a")), # First row with map
  # Second row with histogram and fixed effects
  ggpubr::ggarrange(FigS1b_StatusQuoBiomassStackedDistByManagePlot, FigS1c_allCoeffPlot, ncol = 2, widths = c(1, 3), labels = c("b", ""), vjust = 1),
  nrow = 2, heights = c(1.2, 1) # Label of the map
)

#Plot as a tiff
ggsave(plot = FigS1_MapSizeColorHistGradCoeff,
       filename = paste0(plotDir, "Caldwelletal_RealizedPotentialGains_FigS1_BiomassMapSizeColorHistGradientCoeffPlots.tiff"),
       dpi = 1000, device = "tiff", bg = "white", width = 13, height = 9)

#Plot as a PDF
ggsave(plot = FigS1_MapSizeColorHistGradCoeff,
       filename = paste0(plotDir, "Caldwelletal_RealizedPotentialGains_FigS1_BiomassMapSizeColorHistGradientCoeffPlots.pdf"),
       dpi = 1000, device = "pdf", bg = "white", width = 13, height = 9)


######    Fig S2. Comparison of residuals with distance to the nearest MPA ####
####  Limit the data to fished sites that are within 50km of the nearest no take MPA
fishedLessThan50kmTBL <- predictTBL %>% 
  filter(Management == "Fished" & DistNearMPA_Km <= 50) %>% #436 sites
  select(-contains("Sample"), -contains("SpaMM")) %>% 
  mutate(logBiomassKgHa = log(BiomassKgHa))

#Make predictions of original biomass (not corrected into staus quo) to compare with surveyed data
i = 1
for(i in 1:nrow(bestSpammResTBL)) {
  message("Started model ", i, " of ", nrow(bestSpammResTBL))
  iterSpaMM <- readRDS(file = paste0(modelDir, bestSpammResTBL$ModelFilename[i]))
  
  fishedLessThan50kmTBL <- fishedLessThan50kmTBL %>% 
    mutate(iterPredLogBiomass = unname(predict(object = iterSpaMM,
                                               binding = NA,
                                               type = "link", 
                                               newdata = fishedLessThan50kmTBL,
                                               re.form = as.formula(gsub(pattern = "+Matern(1|Easting+Northing)",
                                                                         replacement = "",
                                                                         x = bestSpammResTBL$Formula[i],
                                                                         fixed = T)))))
  
  colnames(fishedLessThan50kmTBL) <- gsub(pattern = "iterPredLogBiomass", replacement = paste0("SpaMM_", i, "_PredLogBiomass"), x = colnames(fishedLessThan50kmTBL))
}

#Calculate multi-model mean predictions
fishedLessThan50kmTBL <- fishedLessThan50kmTBL %>% 
  mutate(meanPredLogBiomass = rowMeans(fishedLessThan50kmTBL %>% select(starts_with("SpaMM_")))) %>% 
  mutate(meanResidLogBiomass = logBiomassKgHa - meanPredLogBiomass)

#Run models
fishedLessThan50kmTBL.gam <- gam(meanResidLogBiomass ~ s(DistNearMPA_Km), data = fishedLessThan50kmTBL, method = "REML")

#Summaries
summary(fishedLessThan50kmTBL.gam) #p < 0.001; deviance = 0.192%; Rsq = -0.00118

#Model checks
k.check(fishedLessThan50kmTBL.gam) #9
gam.check(fishedLessThan50kmTBL.gam, pch = 19)
appraise(fishedLessThan50kmTBL.gam)
concurvity(fishedLessThan50kmTBL.gam)
draw(fishedLessThan50kmTBL.gam, residuals = T)

#Create new data and plot
data_gamLessThan50km.list <- with(fishedLessThan50kmTBL, list(DistNearMPA_Km = seq(min(DistNearMPA_Km), 50, len = 100)))
newdata_fishedLessThan50kmGam <- emmeans(fishedLessThan50kmTBL.gam, ~DistNearMPA_Km, at = data_gamLessThan50km.list) %>% 
  as.data.frame()

FigS2_FishedResidVsDistMpaLessThan50kmGamPlot <- ggplot(newdata_fishedLessThan50kmGam, aes(x = DistNearMPA_Km, y = emmean)) +
  geom_point(data = fishedLessThan50kmTBL, aes(x = DistNearMPA_Km, y = meanResidLogBiomass), shape = 21) +
  geom_line() + 
  geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL), fill = "black", alpha = 0.3) +
  scale_y_continuous(name = "Multi-model mean residuals\nlog(Biomass)") +
  scale_x_continuous(name = "Distance to nearest MPA (km)") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw() +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = 'black'),
        aspect.ratio = 1/2)

#Plot as a tiff
ggsave(plot = FigS2_FishedResidVsDistMpaLessThan50kmGamPlot,
       filename = paste0(plotDir, "Caldwelletal_RealizedPotentialGains_FigS2_FishedResidVsDistMpaLessThan50kmGamPlot.tiff"),
       dpi = 1000, device = "tiff", bg = "white", width = 8, height = 6)

#Plot as a PDF
ggsave(plot = FigS2_FishedResidVsDistMpaLessThan50kmGamPlot,
       filename = paste0(plotDir, "Caldwelletal_RealizedPotentialGains_FigS2_FishedResidVsDistMpaLessThan50kmGamPlot.pdf"),
       dpi = 1000, device = "pdf", bg = "white", width = 8, height = 6)

######    Fig S3. Test of size and age break (big and old vs. small or new)  ####
######    Create test breaks and  results table 

######    Prepare data 
#### Set up grouped columns 
responseCol <- "BiomassKgHa"
geogCols <- c("Easting", "Northing")
methodCols <- c("Habitat", "Depth", "CensusMethod", "MPAage", "NTZarea")
humanPredCols <- c("Management")
ecoregionCols <- c("ECOREGION", "PROVINCE")
allCols <- c(responseCol, geogCols, methodCols, humanPredCols, ecoregionCols)

#Extract only the high compliance MPAs and the columns of interest
highComplTBL <- predictTBL %>% 
  filter(Management %in% c("UnfishedHighBigOld", "UnfishedHighSmallNew")) %>% 
  select(all_of(allCols))

#create results table
sizeAgeBreakResTBL <- expand.grid(SizeBreaks = c(10,20,30,40),
                                  AgeBreaks = c(10,15,20,25,30)) %>% 
  mutate(RatioBigOldSmallNew = NA,
         BigOldEst = NA,
         BigOldCondSE = NA,
         BigOldConfInt95 = NA,
         ModelR2m = NA,
         ModelR2c = NA)

modelFormula <- as.formula("BiomassKgHa~Habitat+Depth+CensusMethod+(1|PROVINCE/ECOREGION)+Matern(1|Easting+Northing)+Management")

######    For each break, run spaMM model, get difference between "small or new" and "big and new", and add to table 
sizeAgeBreakTestResFilename <- paste0(dataDir, "Caldwelletal_RealizedPotentialGains_TestingAgeSizeBreaksResults.csv")

if(file.exists(sizeAgeBreakTestResFilename)) {
  message("Results from testing size and age breaks already saved")
  sizeAgeBreakResTBL <- read_csv(sizeAgeBreakTestResFilename)
} else {
  i = 1
  for(i in 1:nrow(sizeAgeBreakResTBL)) {
    message("Started row ", i, " of ", nrow(sizeAgeBreakResTBL))
    #Split the dataset into new categories
    iterhighComplTBL <- highComplTBL %>% 
      mutate(Management = ifelse(test = NTZarea > sizeAgeBreakResTBL$SizeBreaks[i] & MPAage > sizeAgeBreakResTBL$AgeBreaks[i],
                                 yes = "UnfishedHighBigOld",
                                 no = "UnfishedHighSmallNew")) %>% 
      mutate(Management = factor(Management, levels = c("UnfishedHighSmallNew", "UnfishedHighBigOld")))
    
    
    sizeAgeBreakResTBL$RatioBigOldSmallNew[i] <- sum(iterhighComplTBL$Management == "UnfishedHighBigOld")/sum(iterhighComplTBL$Management == "UnfishedHighSmallNew")
    
    #Run spaMM model
    newSpammModel <- spaMM::fitme(formula = modelFormula, data = iterhighComplTBL, family = Gamma(log))
    
    #Get estimates from model
    standBetaTableDF <- spaMM::summary.HLfit(newSpammModel)$beta_table %>% 
      as.data.frame() %>% 
      rownames_to_column() %>% 
      rename(Covariate = rowname,
             CondSE = `Cond. SE`) %>% 
      mutate(ConfInt95 = 1.96 * CondSE)
    
    sizeAgeBreakResTBL$BigOldEst[i] <- standBetaTableDF$Estimate[standBetaTableDF$Covariate == "ManagementUnfishedHighBigOld"]
    sizeAgeBreakResTBL$BigOldCondSE[i] <- standBetaTableDF$CondSE[standBetaTableDF$Covariate == "ManagementUnfishedHighBigOld"]
    sizeAgeBreakResTBL$BigOldConfInt95[i] <- standBetaTableDF$ConfInt95[standBetaTableDF$Covariate == "ManagementUnfishedHighBigOld"]
    
    #Get R2 values
    designMatrix <- get_matrix(newSpammModel)
    fixedEff <- fixef(newSpammModel)
    varFix <- var(as.vector(fixedEff %*% t(designMatrix)))
    varRand <- sum(VarCorr(newSpammModel)[VarCorr(newSpammModel)$Group != "Residual","Variance"])    
    varResid <- VarCorr(newSpammModel)[VarCorr(newSpammModel)$Group == "Residual","Variance"]
    
    sizeAgeBreakResTBL$ModelR2m[i] <- varFix/(varFix + varRand + varResid)
    sizeAgeBreakResTBL$ModelR2c[i] <- (varFix + varRand)/(varFix + varRand + varResid)
    
  }
  
  ######    Save overall results, including the R^2 value 
  write_csv(x = sizeAgeBreakResTBL, file = sizeAgeBreakTestResFilename)
}

######    Plot the effect size (difference in biomass) across size and age breaks to identify the best ones
#Add labels to each of the rows of the dataframe
sizeAgeBreakResTBL <- sizeAgeBreakResTBL %>% 
  arrange(RatioBigOldSmallNew) %>% 
  mutate(Label = paste0(SizeBreaks, " km^2 & ", AgeBreaks, " years old (", round(RatioBigOldSmallNew, digits = 2), ")")) %>% 
  mutate(Label = factor(Label, levels = Label))

FigS3_SizeAgeBreaksPlot <- ggplot(sizeAgeBreakResTBL,
                                  mapping = aes(x = Label,
                                                y = BigOldEst,
                                                ymin = BigOldEst - BigOldCondSE,
                                                ymax = BigOldEst + BigOldCondSE)) +
  geom_pointrange() +
  geom_point(aes(x = 20, y = BigOldEst[SizeBreaks == 10 & AgeBreaks == 10]),
             fill = "yellow", shape = 21, size = 3) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_hline(yintercept = 0, color = "dark gray") +
  geom_vline(xintercept = 15.5, linetype = "dashed", color = "gray") +
  xlab("") +
  ylab("Difference in log fish biomass\n(Big and Old vs. Small or New)") +
  coord_flip()

#Plot as a tiff
ggsave(plot = FigS3_SizeAgeBreaksPlot,
       filename = paste0(plotDir, "Caldwelletal_RealizedPotentialGains_FigS3_ManagementSizeAgeBreakPlot.tiff"),
       dpi = 1000, device = "tiff", bg = "white", width = 6, height = 6)

#Plot as a PDF
ggsave(plot = FigS3_SizeAgeBreaksPlot,
       filename = paste0(plotDir, "Caldwelletal_RealizedPotentialGains_FigS3_ManagementSizeAgeBreakPlot.pdf"),
       dpi = 1000, device = "pdf", bg = "white", width = 6, height = 6)


######    Fig S4. Random effects in the most predictive spatial GLMMs  ####
######        S4a. Spatial autocorrelation plot (could make this for all the top models) (based on https://datascienceplus.com/spatial-regression-in-r-part-1-spamm-vs-glmmtmb/) #####
#Get the autocorrelation parameters from each of the best models
spatAutoParamTBL <- data.frame(ModelNum = as.numeric(NA),
                               nu = as.numeric(NA),
                               rho = as.numeric(NA)) %>% 
  filter(!is.na(ModelNum))

spatAutoValsTBL <- data.frame(ModelNum = as.numeric(NA),
                              DistPairs = as.numeric(NA),
                              CorrPairs = as.numeric(NA)) %>% 
  filter(!is.na(ModelNum))

i = 1
for(i in 1:nrow(bestSpammResTBL)) {
  message("Started model ", i)
  #Open the best model
  iterBestStandSpammModel <- readRDS(file = paste0(modelDir, bestSpammResTBL$ModelFilename[i]))
  
  iterNu <- iterBestStandSpammModel$ranFix$corrPars$`3`$nu
  iterRho <- iterBestStandSpammModel$ranFix$corrPars$`3`$rho
  
  spatAutoParamTBL[i,] <- c(i, iterNu, iterRho)
  
  #Create a distance matrix of all sites
  dd <- dist(predictTBL[,c("Easting", "Northing")])
  mm <- MaternCorr(dd, nu = iterNu, rho = iterRho)
  
  iterSpatAutoValsTBL <- data.frame(DistPairs = round(as.numeric(dd)/1000),
                                    CorrPairs = round(as.numeric(mm), digits = 3),
                                    ModelNum = i) %>% 
    filter(DistPairs > 0 & CorrPairs > 0.01) %>% 
    arrange(DistPairs) %>% 
    distinct()
  
  spatAutoValsTBL <- rbind(spatAutoValsTBL, iterSpatAutoValsTBL)
}

nuRangeChar <- paste0(signif(min(spatAutoParamTBL$nu), digits = 3), " - ", signif(max(spatAutoParamTBL$nu), digits = 3))
rhoRangeChar <- paste0(signif(min(spatAutoParamTBL$rho), digits = 3), " -\n           ", signif(max(spatAutoParamTBL$rho), digits = 3))

FigS4a_SpatAutoPlot <- ggplot(data = spatAutoValsTBL, aes(x = DistPairs, y = CorrPairs, group = ModelNum)) +
  geom_line() +
  scale_x_continuous(name = "Distance between pairs of sites (km)",
                     expand = c(0,0),
                     limits = c(0,NA)) +
  scale_y_continuous(name = "Estimated correlation",
                     expand = c(0,0),
                     limits = c(0,NA),
                     breaks = c(0,0.2,0.4,0.6,0.8)) +
  geom_text(x = 100, y = 0.8, label = paste0("nu = ", nuRangeChar,
                                             "\nrho = ", rhoRangeChar), hjust = 0) +
  theme(plot.background = element_blank() ,
        panel.grid.major = element_blank() ,
        panel.grid.minor = element_blank() ,
        panel.border = element_blank() ,
        panel.background = element_blank(),
        axis.line = element_line(color = 'black'),
        aspect.ratio = 1)

######        S4b. Plot showing all random effects of Province and Ecoregion within Province ####
####Create a dataset of the predictions with uncertainties
###Ecoregion
newEcoregionDataTBL <- predictTBL %>% 
  dplyr::select(REALM, PROVINCE, ECOREGION, Long_DD, Lat_DD) %>% 
  arrange(Long_DD, desc(Lat_DD)) %>% 
  dplyr::select(REALM, PROVINCE, ECOREGION) %>% 
  distinct() %>% 
  mutate(REALM = factor(REALM, levels = unique(REALM)),
         PROVINCE = factor(PROVINCE, levels = unique(PROVINCE)),
         Habitat = "Slope",
         Depth = "4-10m",
         CensusMethod = "Belt transect",
         Management = "Fished",
         logGrav_NearMarket = mean(predictTBL$logGrav_NearMarket),
         logGrav_NearPop = mean(predictTBL$logGrav_NearPop),
         logDHW_max_2yr = mean(predictTBL$logDHW_max_2yr),
         logChlA_mean_2yr = mean(predictTBL$logChlA_mean_2yr),
         logChlA_max_2yr = mean(predictTBL$logChlA_max_2yr),
         PAR_mean_2yr = mean(predictTBL$PAR_mean_2yr),
         PAR_sd_2yr = mean(predictTBL$PAR_sd_2yr),
         PAR_skewness_2yr = mean(predictTBL$PAR_skewness_2yr),
         logwave_mean = mean(predictTBL$logwave_mean),
         SST_mean_2yr = mean(predictTBL$SST_mean_2yr),
         SST_max_2yr = mean(predictTBL$SST_max_2yr),
         SST_skewness_2yr = mean(predictTBL$SST_skewness_2yr),
         logPAR_kurtosis_2yr = mean(predictTBL$logPAR_kurtosis_2yr)) 

i = 1
for(i in 1:nrow(bestSpammResTBL)) {
  message("Started model ", i)
  #Open the best model
  iterBestStandSpammModel <- readRDS(file = paste0(modelDir, bestSpammResTBL$ModelFilename[i]))
  
  iterNewEcoregionDataDF <- as.data.frame(get_intervals(iterBestStandSpammModel, 
                                                        newdata = newEcoregionDataTBL, 
                                                        intervals = "predVar",
                                                        type = "link",
                                                        re.form = ~ 1 + (1 | PROVINCE/ECOREGION))) %>% 
    mutate(PredLogBiomass = unname(predict(object = iterBestStandSpammModel,
                                           type = "link", #Predict on log scale for later calculating averages
                                           binding = NA,
                                           newdata = newEcoregionDataTBL,
                                           re.form = ~ 1 + (1 | PROVINCE/ECOREGION))),
           EffectSizeDiff025 = PredLogBiomass - predVar_0.025,
           EffectSizeDiff975 = predVar_0.975 - PredLogBiomass,
           SqEffectSizeDiff025 = EffectSizeDiff025^2,
           SqEffectSizeDiff975 = EffectSizeDiff975^2)
  
  colnames(iterNewEcoregionDataDF) <- paste0("Model", i, "_", colnames(iterNewEcoregionDataDF))
  
  newEcoregionDataTBL <- bind_cols(newEcoregionDataTBL, iterNewEcoregionDataDF)
}

#### Also need to add the among model variance (95% CIs) ####
newEcoregionDataTBL <- newEcoregionDataTBL %>%
  rowwise() %>%
  mutate(AllModels_WithinModelsPropEffectSizeDiff025 = sqrt(mean(as.numeric(across(contains("SqEffectSizeDiff025"))), na.rm = T)),
         AllModels_WithinModelsPropEffectSizeDiff975 = sqrt(mean(as.numeric(across(contains("SqEffectSizeDiff975"))), na.rm = T))) %>% 
  ungroup() %>%
  mutate(ECOREGION_PROVINCE = paste0(as.character(ECOREGION), ":", as.character(PROVINCE)))

#Get the random effect sizes from the model
i = 1
for(i in 1:nrow(bestSpammResTBL)) {
  message("Started model ", i)
  #Open the best model
  iterBestStandSpammModel <- readRDS(file = paste0(modelDir, bestSpammResTBL$ModelFilename[i]))
  
  iterRanef <- ranef(iterBestStandSpammModel)
  
  iterEcoregionRanefDF <- data.frame(EffectSize = unname(iterRanef$`( 1 | ECOREGION:PROVINCE )`),
                                     ECOREGION_PROVINCE = names(iterRanef$`( 1 | ECOREGION:PROVINCE )`))
  
  colnames(iterEcoregionRanefDF)[colnames(iterEcoregionRanefDF) == "EffectSize"] <- paste0("Model", i, "_EffectSize")
  
  newEcoregionDataTBL <- newEcoregionDataTBL %>% 
    left_join(iterEcoregionRanefDF, by = "ECOREGION_PROVINCE") 
}

newEcoregionDataTBL <- newEcoregionDataTBL %>%
  rowwise() %>%
  mutate(AllModelsRandEffectSize = mean(as.numeric(across(ends_with("_EffectSize"))), na.rm = T),
         AmongModelsRandEffectConfInt95 = std.error(as.numeric(across(contains("_EffectSize"))))*1.96) %>% 
  ungroup() %>% 
  mutate(AllModelsRandTotalPropEffectSize025 = sqrt((AllModels_WithinModelsPropEffectSizeDiff025^2) + (AmongModelsRandEffectConfInt95^2)),
         AllModelsRandTotalPropEffectSize975 = sqrt((AllModels_WithinModelsPropEffectSizeDiff975^2) + (AmongModelsRandEffectConfInt95^2)),
         AllModelsRandEffectSize025 = AllModelsRandEffectSize - AllModelsRandTotalPropEffectSize025,
         AllModelsRandEffectSize975 = AllModelsRandEffectSize + AllModelsRandTotalPropEffectSize975)

##Plot the ecoregion effects by province
FigS4b_EcoregionStdEffPlot <- ggplot(newEcoregionDataTBL, 
                                     mapping = aes(x = ECOREGION, 
                                                   y = AllModelsRandEffectSize,
                                                   ymin = AllModelsRandEffectSize025,
                                                   ymax = AllModelsRandEffectSize975)) +
  facet_grid(rows = vars(PROVINCE), scales = "free", space = "free") +
  geom_pointrange(colour = "gray") +
  geom_hline(yintercept = 0, color = "dark gray") +
  scale_y_continuous(limits = c(-1,1)) +
  xlab("Marine Ecoregion") +
  ylab("Standardized Effect Size") +
  labs(tag = "Marine Province") +
  coord_flip() +
  theme(#aspect.ratio = 0.1,
    legend.position = "none",
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    strip.text.y.right = element_text(angle = 0),
    panel.spacing = unit(0.1, "lines"),
    plot.tag = element_text(angle = -90),
    plot.tag.position = c(1.02, 0.5),
    plot.margin = unit(c(1,2,1,1), "lines"))


FigS4_RandomEffects <- ggarrange(FigS4a_SpatAutoPlot,
                                 FigS4b_EcoregionStdEffPlot,
                                 nrow = 2,
                                 labels = c("a", "b"),
                                 heights = c(1,3),
                                 widths = c(1,1))

#Plot as a tiff
ggsave(plot = FigS4_RandomEffects,
       filename = paste0(plotDir, "Caldwelletal_RealizedPotentialGains_FigS4_RandomEffectsPlot.tiff"),
       dpi = 1000, device = "tiff", bg = "white", width = 8, height = 10)

#Plot as a PDF
ggsave(plot = FigS4_RandomEffects,
       filename = paste0(plotDir, "Caldwelletal_RealizedPotentialGains_FigS4_RandomEffectsPlot.pdf"),
       dpi = 1000, device = "pdf", bg = "white", width = 8, height = 10)


######    Fig S5 - Model diagnostics ####
dharmaResidAllModelsTBL <- tibble(Model = as.integer(NA),
                                  DHARMaResiduals = as.numeric(NA)) %>% 
  filter(!is.na(Model))

#Create simulations for for each of the model objects
for(i in 1:nrow(bestSpammResTBL)) {
  #Load the model
  iterSpaMM <- readRDS(file = paste0(modelDir, bestSpammResTBL$ModelFilename[i]))
  
  #Get the DHARMa residual simulations, running how we predicted the biomasses (without re-running the Matern function)
  iterSims <- DHARMa::simulateResiduals(iterSpaMM, re.form = as.formula(paste0("BiomassKgHa~",
                                                                               gsub(pattern = " + Matern(1 | Easting + Northing)",
                                                                                    replacement = "",
                                                                                    fixed = T,
                                                                                    x = as.character(formula(iterSpaMM))[3]))))
  
  #Assign the simulations to the model num
  assign(paste0("spaMM_", i, "_DHARMaSims"), iterSims)
  
  #Extract the residuals and add to data
  dharmaResidAllModelsTBL <- dharmaResidAllModelsTBL %>% 
    bind_rows(tibble(DHARMaResiduals = residuals(iterSims),
                     Model = i))
  
}

### Create a qqplot testing like in DHARMa package but with multi-model means (based on the plotQQunif)
#Calculate multi-model residuals
multiModelScaledResiduals <- rowMeans(tibble(SpaMM_1_ScaledResiduals = spaMM_1_DHARMaSims$scaledResiduals,
                                             SpaMM_2_ScaledResiduals = spaMM_2_DHARMaSims$scaledResiduals,
                                             SpaMM_3_ScaledResiduals = spaMM_3_DHARMaSims$scaledResiduals,
                                             SpaMM_4_ScaledResiduals = spaMM_4_DHARMaSims$scaledResiduals))

#Calculate the mean values for the simulated responses from the DHARMa object
multiModelSimResponse <- (spaMM_1_DHARMaSims$simulatedResponse +
                            spaMM_2_DHARMaSims$simulatedResponse +
                            spaMM_3_DHARMaSims$simulatedResponse +
                            spaMM_4_DHARMaSims$simulatedResponse)/4

#Apply those back to a DHARMa object
multiModelMean_DHARMaSims <- spaMM_1_DHARMaSims
multiModelMean_DHARMaSims$scaledResiduals <- multiModelScaledResiduals
multiModelMean_DHARMaSims$simulatedResponse <- multiModelSimResponse

tiff(filename = paste0(plotDir, "Caldwelletal_RealizedPotentialGains_FigS5_MultiModelDiagnosticPlot.tiff"),
     units = "in",
     width = 8,
     height = 6,
     bg = "white",
     res = 1000)
par(mfrow = c(1,2))
plotQQunif(simulationOutput = multiModelMean_DHARMaSims, testUniformity = F, testOutliers = F, testDispersion = F)
plotResiduals(simulationOutput = multiModelMean_DHARMaSims)
dev.off()

######    Figs S6 & S7. Representativeness of social, ecological, and environmental variables for predictors in the best models  ####
##Subset the reef mask to only include tropical locations and those within the same ecoregions
tropicalLat = 23.43695
reefMaskTBL <- reefMaskTBL %>% 
  filter(Lat_DD < tropicalLat & Lat_DD > -tropicalLat & ECOREGION %in% predictTBL$ECOREGION)  #49878 sites

envPreds <- c("Grav_NearMarket", "Grav_NearPop", 
              "SST_mean_2yr", "SST_max_2yr", "SST_sd_2yr", "SST_kurtosis_2yr", "SST_skewness_2yr", "SSTa_mean_2yr",
              "DHW_mean_2yr", "DHW_max_2yr", "PAR_mean_2yr", "PAR_sd_2yr", "PAR_skewness_2yr", "PAR_kurtosis_2yr",
              "ChlA_mean_2yr", "ChlA_max_2yr", "wave_mean")

envPredLabels <- c(bquote("Nearest market gravity"),
                   bquote("Nearest population gravity"),
                   bquote("Mean SST\n(°C) [1]"),
                   bquote("Maximum SST\n(°C)"),
                   bquote("SST standard deviation"),
                   bquote("SST kurtosis"),
                   bquote("SST skewness"),
                   bquote("Mean SST anomaly\n(°C)"),
                   bquote("Mean DHW\n(°C-weeks) [3]"),
                   bquote("Maximum DHW\n(°C-weeks)"),
                   bquote('Mean PAR\n(Einstein/m^2/day)'),
                   bquote("PAR standard deviation"),
                   bquote("PAR skewness [3]"),
                   bquote("PAR kurtosis"),
                   bquote('Mean Chl-a\n(mg/m^3)'),
                   bquote('Maximum Chl-a\n(mg/m^3)'),
                   bquote("Mean wave energy\n(kW/m) [1]"))

names(envPredLabels) <- envPreds

logEnvPreds <- c("Grav_NearMarket", "Grav_NearPop","SST_kurtosis_2yr", "DHW_mean_2yr", "DHW_max_2yr",
                 "PAR_kurtosis_2yr", "ChlA_mean_2yr", "ChlA_max_2yr", "wave_mean")

#Combine the reef mask predictors with the predict ones in a single dataset
reefMaskEnvPredsTBL <- reefMaskTBL %>% 
  dplyr::select(all_of(envPreds)) %>% 
  mutate(Source = "NOAA Reef Mask Sites")

surveyEnvPredsTBL <- predictTBL %>% 
  mutate(Grav_NearMarket = exp(logGrav_NearMarket),
         Grav_NearPop = exp(logGrav_NearPop),
         SST_kurtosis_2yr = exp(logSST_kurtosis_2yr),
         DHW_mean_2yr = exp(logDHW_mean_2yr),
         DHW_max_2yr = exp(logDHW_max_2yr),
         PAR_kurtosis_2yr = exp(logPAR_kurtosis_2yr),
         ChlA_mean_2yr = exp(logChlA_mean_2yr),
         ChlA_max_2yr = exp(logChlA_max_2yr),
         wave_mean = exp(logwave_mean)) %>% 
  dplyr::select(all_of(envPreds)) %>% 
  mutate(Source = "Survey Sites")

allEnvPredsTBL <- bind_rows(reefMaskEnvPredsTBL, surveyEnvPredsTBL)

#Change all zeros to 0.1 for any of the log transformed variables
allEnvPredsTBL[,logEnvPreds][allEnvPredsTBL[,logEnvPreds] == 0] <- 0.1

#Make a dummy plot to extract the legend
dummyPlot <- ggplot(data = allEnvPredsTBL, aes(x = wave_mean, fill = Source, linetype = Source)) + 
  geom_density(size = 1, alpha = 0.5) +
  scale_x_continuous(name = envPredLabels,
                     expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0),
                     name = "Density") +
  theme(legend.position = "right",
        legend.title = element_blank(),
        plot.background = element_blank() ,
        panel.grid.major = element_blank() ,
        panel.grid.minor = element_blank() ,
        panel.border = element_blank() ,
        panel.background = element_blank(),
        axis.line = element_line(color = 'black'),
        aspect.ratio = 1)

figS6_legend <- get_legend(dummyPlot)

###Create plots for all of the predictors
i = 1
for (i in 1:length(envPreds)) {
  message("Started plotting predictor ", i, " of ", length(envPreds), ": ", envPreds[i])
  envPred = envPreds[i]
  envPredLabel = envPredLabels[[i]]
  
  labelPosition = "none"
  
  allEnvPredsTBL[,"EnvPred"] <- allEnvPredsTBL[,envPred]
  
  if(envPred %in% logEnvPreds) {
    densityPlot <- ggplot(data = allEnvPredsTBL, aes(x = EnvPred, fill = Source, linetype = Source)) + 
      geom_density(size = 1, alpha = 0.5) +
      scale_x_log10(name = envPredLabel,
                    labels = label_number(drop0trailing = T),
                    expand = c(0,0)) +
      scale_y_continuous(expand = c(0,0),
                         name = "Density") +
      theme(legend.position = labelPosition,
            legend.title = element_blank(),
            plot.background = element_blank() ,
            panel.grid.major = element_blank() ,
            panel.grid.minor = element_blank() ,
            panel.border = element_blank() ,
            panel.background = element_blank(),
            axis.line = element_line(color = 'black'),
            aspect.ratio = 1)
  } else {
    densityPlot <- ggplot(data = allEnvPredsTBL, aes(x = EnvPred, fill = Source, linetype = Source)) + 
      geom_density(size = 1, alpha = 0.5) +
      scale_x_continuous(name = envPredLabel,
                         expand = c(0,0)) +
      scale_y_continuous(expand = c(0,0),
                         name = "Density") +
      theme(legend.position = labelPosition,
            legend.title = element_blank(),
            plot.background = element_blank() ,
            panel.grid.major = element_blank() ,
            panel.grid.minor = element_blank() ,
            panel.border = element_blank() ,
            panel.background = element_blank(),
            axis.line = element_line(color = 'black'),
            aspect.ratio = 1)
  }
  
  assign(paste0(envPred, "_ReefMaskVsSurveyDensityPlot"), densityPlot)
}

bestModelPreds <- c("Grav_NearMarket", "Grav_NearPop", 
                    "SST_mean_2yr", "SST_max_2yr", "SST_skewness_2yr", 
                    "DHW_max_2yr", "PAR_mean_2yr", "PAR_sd_2yr", "PAR_skewness_2yr", "PAR_kurtosis_2yr",
                    "ChlA_mean_2yr", "ChlA_max_2yr", "wave_mean")

###Assemble a figure with those predictors that were in the best models
FigS6_RepresentativenessBestModelPredPlots <- ggpubr::ggarrange(Grav_NearMarket_ReefMaskVsSurveyDensityPlot,
                                                                Grav_NearPop_ReefMaskVsSurveyDensityPlot,
                                                                SST_mean_2yr_ReefMaskVsSurveyDensityPlot,
                                                                SST_max_2yr_ReefMaskVsSurveyDensityPlot,
                                                                SST_skewness_2yr_ReefMaskVsSurveyDensityPlot,
                                                                DHW_max_2yr_ReefMaskVsSurveyDensityPlot,
                                                                PAR_mean_2yr_ReefMaskVsSurveyDensityPlot,
                                                                PAR_sd_2yr_ReefMaskVsSurveyDensityPlot,
                                                                PAR_skewness_2yr_ReefMaskVsSurveyDensityPlot,
                                                                PAR_kurtosis_2yr_ReefMaskVsSurveyDensityPlot,
                                                                ChlA_mean_2yr_ReefMaskVsSurveyDensityPlot,
                                                                ChlA_max_2yr_ReefMaskVsSurveyDensityPlot,
                                                                wave_mean_ReefMaskVsSurveyDensityPlot, ncol = 4, nrow = 4, common.legend = T, legend = "bottom")
                                                                

#Plot as a tiff
ggsave(plot = FigS6_RepresentativenessBestModelPredPlots,
       filename = paste0(plotDir, "Caldwelletal_RealizedPotentialGains_FigS6_RepresentativenessBestModelPredPlots.tiff"),
       dpi = 1000, device = "tiff", bg = "white", width = 8.5, height = 12)

#Plot as a PDF
ggsave(plot = FigS6_RepresentativenessBestModelPredPlots,
       filename = paste0(plotDir, "Caldwelletal_RealizedPotentialGains_FigS6_RepresentativenessBestModelPredPlots.pdf"),
       dpi = 1000, device = "pdf", bg = "white", width = 8.5, height = 12)

otherModelPreds <- setdiff(envPreds, bestModelPreds)

FigS7_RepresentativenessOtherModelPredPlots <- ggpubr::ggarrange(SST_sd_2yr_ReefMaskVsSurveyDensityPlot,
                                                                 SST_kurtosis_2yr_ReefMaskVsSurveyDensityPlot,
                                                                 SSTa_mean_2yr_ReefMaskVsSurveyDensityPlot,
                                                                 DHW_mean_2yr_ReefMaskVsSurveyDensityPlot,
                                                                 ncol = 2, nrow = 2, common.legend = T, legend = "bottom")

#Plot as a tiff
ggsave(plot = FigS7_RepresentativenessOtherModelPredPlots,
       filename = paste0(plotDir, "Caldwelletal_RealizedPotentialGains_FigS7_RepresentativenessOtherModelPredPlots.tiff"),
       dpi = 1000, device = "tiff", bg = "white", width = 4, height = 6)

#Plot as a PDF
ggsave(plot = FigS7_RepresentativenessOtherModelPredPlots,
       filename = paste0(plotDir, "Caldwelletal_RealizedPotentialGains_FigS7_RepresentativenessOtherModelPredPlots.pdf"),
       dpi = 1000, device = "pdf", bg = "white", width = 4, height = 6)


