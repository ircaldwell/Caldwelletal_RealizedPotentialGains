###### Caldwell et al. - Main figures for "Protection efforts have resulted in ~10% of existing fish biomass on global coral reefs ####
######  Author: Iain R. Caldwell
######  Date last revised: May 18, 2023
######  Main figures for paper:
######    Fig 1. Realized gains
######      a. Cumulative biomass vs. % of sites for status quo and fished scenario
######      b. Realized gains given subsampled full protection
######      c. Map of realized gains
######      Assemble figure 1
######    Fig 2. Variation in realized gains within protection categories
######      a. Stacked distribution of realized gains (colored by management)
######      b. Realized fish biomass gains along a gradient of nearest market gravity
######      Assemble figure 2
######    Fig 3. Potential gains
######      a. Cumulative potential and realized gains (% of status quo - potential gains positive, realized gains negative)
######      e. Map of potential full protection gains  
######      Assemble figure 3

rm(list = ls()) #remove past stored objects
options(scipen=999) #disable scientific notation

#Load libraries
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(lemon) #for facet_rep_grid function
library(patchwork)
library(ggnewscale)
library(sp) #for SpatialPoints function
library(ggridges)
library(ggh4x)
library(scales) #for removing trailing zeroes in log plots and custom transformations
library(spaMM)

#Set random seed for reproducibility
set.seed(1234)

#Set the number of samples
numSamples = 1000

####Set parameters and directories ####
dateNum = as.character(Sys.Date())
setwd('..')
resultsDir <- paste0(getwd(), "/Caldwelletal_Results/")
plotDir <- paste0(getwd(), "/Caldwelletal_Plots/")
modelDir <- paste0(getwd(), "/Caldwelletal_Models/")

######  Load predictions ####
predictTBL <- readRDS(file = paste0(resultsDir, "Caldwelletal_RealizedPotentialGains_PredictionsFishBiomass_2023-05-21.rds")) 

#Get numbers for each protection type (Table 1)
manageSamples <- as.data.frame(table(predictTBL$Management))

######  Load list with best models ####
bestSpammResTBL <- read_csv(paste0(resultsDir, "Caldwelletal_RealizedPotentialGains_BestSpammSummResWithFilenames_2023-05-21.csv")) 

#get the number of best models
numModels <- nrow(bestSpammResTBL)

######    Fig 1. Realized gains ####
######      > a. Cumulative biomass vs. % of sites for status quo and fished scenario ####
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

###Sort data for figure
cumBiomassTBL <- predictTBL %>% 
  mutate(Management = factor(Management, levels = c("UnfishedHighBigOld", "UnfishedHighSmallNew", "UnfishedLow", "Restricted", "Fished"))) %>% 
  arrange(desc(Management), desc(StatusQuoBiomassMed)) %>% 
  mutate(PercentageSites = c(1:nrow(predictTBL))/nrow(predictTBL)*100) 

cumStatusQuoBiomassSamplesTBL <- data.frame(matrix(nrow = nrow(cumBiomassTBL), ncol = numSamples*numModels))
cumFishedBiomassSamplesTBL <- data.frame(matrix(nrow = nrow(cumBiomassTBL), ncol = numSamples*numModels))

#Calculate the 95% quantiles and quartiles for the cumulative sums from the sample distributions
modelNum = 1
iterNum <- 0
for(modelNum in 1:numModels) {
  message("Started model ", modelNum)
  i = 1
  for(i in 1:numSamples) {
    statusQuoColname <- paste0("SpaMM_", modelNum, "_StatusQuoPredSample_",i)
    fishedColname <- paste0("SpaMM_", modelNum, "_FishedPredSample_",i)
    
    iterNum <- iterNum + 1
    cumStatusQuoBiomassSamplesTBL[,iterNum] <- (cumsum(cumBiomassTBL[,statusQuoColname])/sum(cumBiomassTBL[,statusQuoColname])*100)[,1]
    cumFishedBiomassSamplesTBL[,iterNum] <- (cumsum(cumBiomassTBL[,fishedColname])/sum(cumBiomassTBL[,statusQuoColname])*100)[,1]
    message("Fished percentage for iter ", iterNum, " is ", max(cumFishedBiomassSamplesTBL[,iterNum]))
  }
}

cumStatusQuoBiomassQuantilesTBL <- data.frame(matrix(nrow = nrow(cumBiomassTBL), ncol = 5))
cumFishedBiomassQuantilesTBL <- data.frame(matrix(nrow = nrow(cumBiomassTBL), ncol = 5))

i = 1
for(i in 1:nrow(cumBiomassTBL)) {
  cumStatusQuoBiomassQuantilesTBL[i,] <- unname(quantile(x = cumStatusQuoBiomassSamplesTBL[i,], probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))
  cumFishedBiomassQuantilesTBL[i,] <- unname(quantile(x = cumFishedBiomassSamplesTBL[i,], probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))
}
  
colnames(cumStatusQuoBiomassQuantilesTBL) <- c("CumStatusQuoBiomassCompStatusQuo_025Quant",
                                               "CumStatusQuoBiomassCompStatusQuo_25Quant",
                                               "CumStatusQuoBiomassCompStatusQuo_50Quant",
                                               "CumStatusQuoBiomassCompStatusQuo_75Quant",
                                               "CumStatusQuoBiomassCompStatusQuo_975Quant")

colnames(cumFishedBiomassQuantilesTBL) <- c("CumFishedBiomassCompStatusQuo_025Quant",
                                            "CumFishedBiomassCompStatusQuo_25Quant",
                                            "CumFishedBiomassCompStatusQuo_50Quant",
                                            "CumFishedBiomassCompStatusQuo_75Quant",
                                            "CumFishedBiomassCompStatusQuo_975Quant")

cumBiomassTBL <- cumBiomassTBL %>% bind_cols(cumStatusQuoBiomassQuantilesTBL, cumFishedBiomassQuantilesTBL)

###>> Fig 1a plot ####
Fig1a_CumBiomassVsPerSitesFishedStatusQuoRealGainsPlot <- ggplot(data = cumBiomassTBL, aes(x = PercentageSites)) +
  geom_ribbon(aes(ymin = max(CumFishedBiomassCompStatusQuo_025Quant), ymax = max(CumFishedBiomassCompStatusQuo_975Quant)),
              fill = gainColors["RealizedGains"]) +
  geom_segment(
    aes(y = max(CumFishedBiomassCompStatusQuo_50Quant), yend = max(CumFishedBiomassCompStatusQuo_50Quant)), x = 0, xend = 100, 
    size = 1, 
    linetype = "dashed",
    colour = "black" 
  ) +
  geom_segment(
    y = 100, yend = 100, x = 0, xend = 100, 
    size = 1, 
    linetype = "dashed",
    colour = "black" 
  ) +
  geom_ribbon(data = cumBiomassTBL %>% filter(Management != "Fished"),
              aes(ymin = CumFishedBiomassCompStatusQuo_025Quant,
                  ymax = CumFishedBiomassCompStatusQuo_975Quant),
              fill = manageColors["Fished"], alpha = 0.5) +
  geom_ribbon(aes(ymin = CumStatusQuoBiomassCompStatusQuo_025Quant,
                  ymax = CumStatusQuoBiomassCompStatusQuo_975Quant,
                  fill = Management),
              alpha = 0.5) +
  geom_line(data = cumBiomassTBL %>% filter(Management != "Fished"),
            aes(y = CumFishedBiomassCompStatusQuo_50Quant),
            color = manageColors["Fished"], linewidth = 1, linetype = "dashed") +
  geom_line(aes(y = CumStatusQuoBiomassCompStatusQuo_50Quant, color = Management), linewidth = 1) +
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
  annotate(geom = "text", x = 79, y = 90, label = "total\nrealized\ngains", color = "black", hjust = 1, vjust = 0.5) +
  #Add arrow for total realized gains
  geom_segment(
    aes(x = 80, y = max(CumFishedBiomassCompStatusQuo_50Quant)), xend = 80, yend = 100,
    lineend = "round", 
    linejoin = "round",
    size = 1, 
    arrow = arrow(length = unit(0.3, "cm")),
    colour = "black" 
  ) +
  #Add arrow and text for the status quo biomass
  geom_segment(
    aes(x = 59, y = 67, xend = 64, yend = 62),
    lineend = "round", 
    linejoin = "round",
    size = 1, 
    arrow = arrow(length = unit(0.3, "cm")),
    colour = "black" 
  ) + 
  annotate(geom = "text", x = 50, y = 70, label = "Status quo", color = "black") +
  #Add arrow pointing to the fished scenario
  geom_segment(
    aes(x = 80, y = 50, xend = 75, yend = 59),
    lineend = "round", 
    linejoin = "round",
    size = 1, 
    arrow = arrow(length = unit(0.3, "cm")),
    colour = "black" 
  ) + 
  annotate(geom = "text", x = 85, y = 48, label = "Fished scenario", color = "black")

######      > b. Realized gains given subsampled full protection ####
# Create new dataset for figure with actual data and subsampled data (1000 subsamples) 
#Get the percentages for the actual data
realGainsFullMpaByPercFilename <- paste0(resultsDir, "RealizedGains_SubsampledByPercFullMPA.csv")

if(file.exists(realGainsFullMpaByPercFilename)) {
  message("Realized gains already subsampled")
  realGainsFullMpaByPercProtTBL <- read_csv(realGainsFullMpaByPercFilename)
} else {
  realGainsFullMpaPercStatusQuoActual <- c()
  
  modelNum <- 1 #for testing
  for(modelNum in 1:numModels) {
    message("Started model ", modelNum)
    i = 1
    for(i in 1:numSamples) {
      realGainsFullMpaPercStatusQuoActual <- c(realGainsFullMpaPercStatusQuoActual,
                                               sum(predictTBL[,paste0("SpaMM_", modelNum, "_RealizedGainsSample_",i)])/
                                                 sum(predictTBL[,paste0("SpaMM_", modelNum, "_StatusQuoPredSample_",i)]))
    }
  }
  
  realGainsFullMpaByPercProtTBL <- data.frame(NumFullMPA = sum(!predictTBL$Management %in% c("Fished", "Restricted")),
                                              PercFullMPA = sum(!predictTBL$Management %in% c("Fished", "Restricted"))/nrow(predictTBL)*100,
                                              RealGainsFullMpaPercStatusQuoQuant50 = unname(quantile(x = realGainsFullMpaPercStatusQuoActual, probs = 0.5)),
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
  
  for(modelNum in 1:numModels) {
    message("Started model ", modelNum)
    j = 1
    for(j in 1:numSamples) {
      realGainsFullMpaPercStatusQuoGlobalProt <- c(realGainsFullMpaPercStatusQuoGlobalProt,
                                                   sum(noFullProtTBL[,paste0("SpaMM_", modelNum, "_RealizedGainsSample_",j)])/
                                                     sum(noFullProtTBL[,paste0("SpaMM_", modelNum, "_StatusQuoPredSample_",j)]))
    }
  }
  
  realGainsFullMpaByPercProtTBL$RealGainsFullMpaPercStatusQuoQuant50[nrow(realGainsFullMpaByPercProtTBL)] <- unname(quantile(x = realGainsFullMpaPercStatusQuoGlobalProt, probs = 0.5))
  realGainsFullMpaByPercProtTBL$RealGainsFullMpaPercStatusQuoQuant975[nrow(realGainsFullMpaByPercProtTBL)] <- unname(quantile(x = realGainsFullMpaPercStatusQuoGlobalProt, probs = 0.975))
  realGainsFullMpaByPercProtTBL$RealGainsFullMpaPercStatusQuoQuant025[nrow(realGainsFullMpaByPercProtTBL)] <- unname(quantile(x = realGainsFullMpaPercStatusQuoGlobalProt, probs = 0.025))
  realGainsFullMpaByPercProtTBL$Data[nrow(realGainsFullMpaByPercProtTBL)] <- "No Full MPAs"
  
  
  #Successively subsample the MPAs and recalculate the realized gains
  subsamples <- seq(from = 1, to = sum(predictTBL$Management %in% c("UnfishedHighBigOld", "UnfishedHighSmallNew", "UnfishedLow"))-1, by = 1)
  
  i = 1
  for(i in subsamples) {
    message("Started subsampling ", i, " full MPAs")
    realGainsFullMpaByPercProtTBL[nrow(realGainsFullMpaByPercProtTBL)+1,] <- NA
    realGainsFullMpaByPercProtTBL$NumFullMPA[nrow(realGainsFullMpaByPercProtTBL)] <- i
    realGainsFullMpaByPercProtTBL$PercFullMPA[nrow(realGainsFullMpaByPercProtTBL)] <- i/(nrow(noFullProtTBL)+i)*100
    
    realGainsFullMpaPercStatusQuoGlobalProt <- c()
    
    for(modelNum in 1:numModels) {
      message("Started model ", modelNum)
      j = 1
      for(j in 1:numSamples) {
        fullRowSample <- sample(x = row.names(fullProtTBL), size = i)
        subsampleTBL <- bind_rows(fullProtTBL[fullRowSample,],
                                  noFullProtTBL)
        realGainsFullMpaPercStatusQuoGlobalProt <- c(realGainsFullMpaPercStatusQuoGlobalProt,
                                                     sum(subsampleTBL[,paste0("SpaMM_", modelNum, "_RealizedGainsSample_",j)])/
                                                       sum(subsampleTBL[,paste0("SpaMM_", modelNum, "_StatusQuoPredSample_",j)]))
      }
    }
    
    realGainsFullMpaByPercProtTBL$RealGainsFullMpaPercStatusQuoQuant50[nrow(realGainsFullMpaByPercProtTBL)] <- unname(quantile(x = realGainsFullMpaPercStatusQuoGlobalProt, probs = 0.5))
    realGainsFullMpaByPercProtTBL$RealGainsFullMpaPercStatusQuoQuant975[nrow(realGainsFullMpaByPercProtTBL)] <- unname(quantile(x = realGainsFullMpaPercStatusQuoGlobalProt, probs = 0.975))
    realGainsFullMpaByPercProtTBL$RealGainsFullMpaPercStatusQuoQuant025[nrow(realGainsFullMpaByPercProtTBL)] <- unname(quantile(x = realGainsFullMpaPercStatusQuoGlobalProt, probs = 0.025))
    realGainsFullMpaByPercProtTBL$Data[nrow(realGainsFullMpaByPercProtTBL)] <- "Subsample"
  }
  
  #Save the dataset so this doesn't have to be run again
  write_csv(x = realGainsFullMpaByPercProtTBL, file = realGainsFullMpaByPercFilename)
}

realGainsFullMpaByPercProtTBL <- realGainsFullMpaByPercProtTBL %>% 
  arrange(NumFullMPA)

#Get the closest percentage to the Allen Coral percentages (3.06%)
closestPercAllenCoralProt <- realGainsFullMpaByPercProtTBL$PercFullMPA[which.min(abs(realGainsFullMpaByPercProtTBL$PercFullMPA - 3.06))] 

###>> Fig 1b plot ####
Fig1b_RealGainsSubsampleAllenCoralPercPlot <- ggplot(data = realGainsFullMpaByPercProtTBL,
                                               aes(x = PercFullMPA,
                                                   y = RealGainsFullMpaPercStatusQuoQuant50*100,
                                                   ymin = RealGainsFullMpaPercStatusQuoQuant025*100,
                                                   ymax = RealGainsFullMpaPercStatusQuoQuant975*100)) +
  geom_line(size = 1) +
  geom_ribbon(alpha = 0.5) +
  geom_segment(data = realGainsFullMpaByPercProtTBL %>% filter(Data == "All surveys"),
               aes(x = 0, xend = PercFullMPA, y = RealGainsFullMpaPercStatusQuoQuant50*100, yend = RealGainsFullMpaPercStatusQuoQuant50*100),
               linetype = "dashed") +
  geom_segment(data = realGainsFullMpaByPercProtTBL %>% filter(Data == "All surveys"),
               aes(x = PercFullMPA, xend = PercFullMPA, y = 5, yend = RealGainsFullMpaPercStatusQuoQuant50*100),
               linetype = "dashed") +
  geom_errorbar(data = realGainsFullMpaByPercProtTBL %>% filter(Data == "All surveys")) +
  geom_point(data = realGainsFullMpaByPercProtTBL %>% filter(Data == "All surveys"), size = 3, shape = 21, fill = "grey") +
  geom_segment(data = realGainsFullMpaByPercProtTBL %>% filter(PercFullMPA == closestPercAllenCoralProt),
               aes(x = 0, xend = PercFullMPA, y = RealGainsFullMpaPercStatusQuoQuant50*100, yend = RealGainsFullMpaPercStatusQuoQuant50*100),
               linetype = "dotted") +
  geom_segment(data = realGainsFullMpaByPercProtTBL %>% filter(PercFullMPA == closestPercAllenCoralProt),
               aes(x = PercFullMPA, xend = PercFullMPA, y = 5, yend = RealGainsFullMpaPercStatusQuoQuant50*100),
               linetype = "dotted") +
  geom_errorbar(data = realGainsFullMpaByPercProtTBL %>% filter(PercFullMPA == closestPercAllenCoralProt)) +
  geom_point(data = realGainsFullMpaByPercProtTBL %>% filter(PercFullMPA == closestPercAllenCoralProt), size = 3, shape = 21, fill = "white") +
  scale_y_continuous(name = "Realized gains\n(% of total status quo biomass)",
                     expand = c(0,0)
  ) +
  scale_x_continuous(name = "% fully protected MPA coverage",
                     expand = c(0,0)) +
  theme_classic() +
  theme(aspect.ratio = 1) 

######      > c. Map of realized gains ####
###Set the minimum value for the log gains
zeroLogTrans = 1

###Set up map
mapWorld <- map_data('world2')

###Prepare data
realGainMapTBL <- predictTBL %>%
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

ggsave(plot = Fig1_RealGainsCumBiomassCorrSubsampleMap, 
       filename = paste0(plotDir, "RealizedGainsMS_Fig1_RealGainsCumBiomassCorrSampleMap_", dateNum, ".tiff"),
       width = 10,
       height = 10)


######    Fig 2. Variation in realized gains within protection categories ####
######      > a. Stacked distribution of realized gains (colored by management) ####
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

modelNum <- 1 #for testing
for(modelNum in 1:numModels) {
  message("Started model ", modelNum)
  i = 1
  for(i in 1:numSamples) {
    restRG <- c(restRG, unname(unlist(predictTBL[predictTBL$Management == "Restricted", paste0("SpaMM_", modelNum, "_RealizedGainsSample_", i)])))
    unfishedLowRG <- c(unfishedLowRG, unname(unlist(predictTBL[predictTBL$Management == "UnfishedLow", paste0("SpaMM_", modelNum, "_RealizedGainsSample_", i)])))
    unfishedHighSmallNewRG <- c(unfishedHighSmallNewRG, unname(unlist(predictTBL[predictTBL$Management == "UnfishedHighSmallNew", paste0("SpaMM_", modelNum, "_RealizedGainsSample_", i)])))
    unfishedHighBigOldRG <- c(unfishedHighBigOldRG, unname(unlist(predictTBL[predictTBL$Management == "UnfishedHighBigOld", paste0("SpaMM_", modelNum, "_RealizedGainsSample_", i)])))
    allRG <- c(allRG, unname(unlist(predictTBL[predictTBL$Management != "Fished", paste0("SpaMM_", modelNum, "_RealizedGainsSample_", i)])))
  }
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

sumCounts <- nrow(bestSpammResTBL)*numSamples*nrow(predictTBL)

freqBreaks <- c(0, 0.25*sumCounts, 0.5*sumCounts, 0.75*sumCounts)
freqLabels <- c(0, 0.25, 0.5, 0.75)

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
        axis.ticks.y = element_blank(),
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

######      > b. Realized fish biomass gains along a gradient of nearest market gravity ####
#Get predictions from each of the best models across the full gradient of gravity
gravGradientUnfishedBigOldTBL <- tibble(logGrav_NearMarket = seq(from = min(predictTBL$logGrav_NearMarket),
                                                                 to = max(predictTBL$logGrav_NearMarket),
                                                                 length.out = 1000)) %>% 
  mutate(Habitat = "Slope",
         Depth = "4-10m",
         CensusMethod = "Belt transect",
         Management = "UnfishedHighBigOld",
         SST_mean_2yr = mean(predictTBL$SST_mean_2yr),
         SST_min_2yr = mean(predictTBL$SST_min_2yr),
         SST_max_2yr = mean(predictTBL$SST_max_2yr),
         SST_sd_2yr = mean(predictTBL$SST_sd_2yr),
         SST_skewness_2yr = mean(predictTBL$SST_skewness_2yr),
         SST_kurtosis_2yr = mean(predictTBL$SST_kurtosis_2yr),
         SSTa_mean_2yr = mean(predictTBL$SSTa_mean_2yr),
         SSTa_min_2yr = mean(predictTBL$SSTa_min_2yr),
         SSTa_max_2yr = mean(predictTBL$SSTa_max_2yr),
         SSTa_sd_2yr = mean(predictTBL$SSTa_sd_2yr),
         SSTa_skewness_2yr = mean(predictTBL$SSTa_skewness_2yr),
         SSTa_kurtosis_2yr = mean(predictTBL$SSTa_kurtosis_2yr),
         logDHW_mean_2yr = mean(predictTBL$logDHW_mean_2yr),
         logDHW_max_2yr = mean(predictTBL$logDHW_max_2yr),
         logDHW_kurtosis_2yr = mean(predictTBL$logDHW_kurtosis_2yr),
         logDHW_skewness_2yr = mean(predictTBL$logDHW_skewness_2yr),
         PAR_mean_2yr = mean(predictTBL$PAR_mean_2yr),
         PAR_min_2yr = mean(predictTBL$PAR_min_2yr),
         PAR_max_2yr = mean(predictTBL$PAR_max_2yr),
         PAR_sd_2yr = mean(predictTBL$PAR_sd_2yr),
         PAR_skewness_2yr = mean(predictTBL$PAR_skewness_2yr),
         PAR_kurtosis_2yr = mean(predictTBL$PAR_kurtosis_2yr),
         logChlA_mean_2yr = mean(predictTBL$logChlA_mean_2yr),
         logChlA_min_2yr = mean(predictTBL$logChlA_min_2yr),
         logChlA_max_2yr = mean(predictTBL$logChlA_max_2yr),
         logChlA_sd_2yr = mean(predictTBL$logChlA_sd_2yr),
         logChlA_skewness_2yr = mean(predictTBL$logChlA_skewness_2yr),
         logChlA_kurtosis_2yr = mean(predictTBL$logChlA_kurtosis_2yr),
         logwave_mean = mean(predictTBL$logwave_mean),
         logGrav_NearPop = mean(predictTBL$logGrav_NearPop))

gravGradientFishedTBL <- gravGradientUnfishedBigOldTBL %>% 
  mutate(Management = "Fished")

gravGradientPredsTBL <- tibble(logGrav_NearMarket = as.numeric(NA),
                               ModelNum = as.numeric(NA),
                               UnfishedBigOldPredMean = as.numeric(NA),
                               UnfishedBigOldPredSD = as.numeric(NA),
                               FishedPredMean = as.numeric(NA),
                               FishedPredSD = as.numeric(NA)) %>% 
  filter(!is.na(ModelNum))

i = 1
for(i in 1:nrow(bestSpammResTBL)) {
  message("Started model ", i, " of ", nrow(bestSpammResTBL))
  iterSpaMM <- readRDS(file = paste0(modelDir, bestSpammResTBL$ModelFilename[i]))
  
  iterGravGradientPredsTBL <- gravGradientUnfishedBigOldTBL %>% 
    select(logGrav_NearMarket) %>% 
    distinct() %>% 
    mutate(ModelNum = i,
           UnfishedBigOldPredMean = unname(predict(object = iterSpaMM,
                                                   binding = NA,
                                                   newdata = gravGradientUnfishedBigOldTBL,
                                                   re.form = NA)),
           UnfishedBigOldPredSD = sqrt(unname(get_predVar(object = iterSpaMM,
                                                          newdata = gravGradientUnfishedBigOldTBL, 
                                                          re.form = NA))),
           FishedPredMean = unname(predict(object = iterSpaMM,
                                           binding = NA,
                                           newdata = gravGradientFishedTBL,
                                           re.form = NA)),
           FishedPredSD = sqrt(unname(get_predVar(object = iterSpaMM,
                                                  newdata = gravGradientFishedTBL, 
                                                  re.form = NA))))
  gravGradientPredsTBL <- rbind(gravGradientPredsTBL, iterGravGradientPredsTBL)
}

gravGradientPredsTBL <- gravGradientPredsTBL %>% 
  mutate(RealizedGainsMean = exp(UnfishedBigOldPredMean) - exp(FishedPredMean),
         ModelNum = as.factor(ModelNum))

#Get the average biomasses for each gravity
gravGradientAllModelsTBL <- gravGradientPredsTBL %>% 
  group_by(logGrav_NearMarket) %>% 
  summarise(UnfishedBigOldPredMean = mean(UnfishedBigOldPredMean),
            FishedPredMean = mean(FishedPredMean)) %>% 
  mutate(RealizedGainsMean = exp(UnfishedBigOldPredMean) - exp(FishedPredMean))

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
  nrow = 2, heights = c(2.1,0.1)
)

ggsave(plot = Fig2_RealGainsStackedHistGradientGravity, 
       filename = paste0(plotDir, "RealizedGainsMS_Fig2_RealGainsStackedHistGradientGravity_", dateNum, ".tiff"),
       width = 10,
       height = 6)

######    Fig 3. Potential gains ####
######      > a. Cumulative potential and realized gains (% of status quo - potential gains positive, realized gains negative) ####
###Sort data for figure
cumGainsTBL <- bind_rows(predictTBL %>%
                           dplyr::filter(Management %in% c("UnfishedHighBigOld", "UnfishedHighSmallNew")) %>% 
                           arrange(desc(Management), desc(RealizedGainMed)),
                         predictTBL %>% 
                           dplyr::filter(!Management %in% c("UnfishedHighBigOld", "UnfishedHighSmallNew")) %>% 
                           arrange(desc(PotentialGainMed))) %>% 
  mutate(PercentageSites = c(1:nrow(predictTBL))/nrow(predictTBL)*100) 
  
cumPotentialGainsPercStatusQuoTBL <- data.frame(matrix(nrow = nrow(cumGainsTBL), ncol = numSamples*numModels))
cumRealizedGainsPercStatusQuoTBL <- data.frame(matrix(nrow = nrow(cumGainsTBL), ncol = numSamples*numModels))

#Calculate the 95% quantiles and quartiles for the cumulative sums from the sample distributions
modelNum <- 1
iterNum <- 0
for(modelNum in 1:numModels) {
  message("Started model ", modelNum)
  i <- 1
  for(i in 1:numSamples) {
    potentialGainsColname <- paste0("SpaMM_", modelNum, "_PotentialGainsSample_",i)
    realizedGainsColname <- paste0("SpaMM_", modelNum, "_RealizedGainsSample_",i)
    statusQuoColname <- paste0("SpaMM_", modelNum, "_StatusQuoPredSample_",i)
    
    iterNum <- iterNum + 1
    cumPotentialGainsPercStatusQuoTBL[,iterNum] <- (cumsum(cumGainsTBL[,potentialGainsColname])/sum(cumGainsTBL[,statusQuoColname])*100)[,1]
    cumRealizedGainsPercStatusQuoTBL[,iterNum] <- (cumsum(cumGainsTBL[,realizedGainsColname])/sum(cumGainsTBL[,statusQuoColname])*100)[,1]*-1
    message("Maximum cumulative gains = ", max(cumPotentialGainsPercStatusQuoTBL[,iterNum]), " (P); ", min(cumRealizedGainsPercStatusQuoTBL[,iterNum]), " (R)")
  }
}

cumPotentialGainsPercStatusQuoQuantilesTBL <- data.frame(matrix(nrow = nrow(cumGainsTBL), ncol = 5))
cumRealizedGainsPercStatusQuoQuantilesTBL <- data.frame(matrix(nrow = nrow(cumGainsTBL), ncol = 5))

i = 1
for(i in 1:nrow(cumGainsTBL)) {
  cumPotentialGainsPercStatusQuoQuantilesTBL[i,] <- unname(quantile(x = cumPotentialGainsPercStatusQuoTBL[i,], probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))
  cumRealizedGainsPercStatusQuoQuantilesTBL[i,] <- unname(quantile(x = cumRealizedGainsPercStatusQuoTBL[i,], probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))
}

colnames(cumPotentialGainsPercStatusQuoQuantilesTBL) <- c("CumPotentialGainsPercStatusQuo_025Quant",
                                                          "CumPotentialGainsPercStatusQuo_25Quant",
                                                          "CumPotentialGainsPercStatusQuo_50Quant",
                                                          "CumPotentialGainsPercStatusQuo_75Quant",
                                                          "CumPotentialGainsPercStatusQuo_975Quant")

colnames(cumRealizedGainsPercStatusQuoQuantilesTBL) <- c("CumRealizedGainsPercStatusQuo_025Quant",
                                                         "CumRealizedGainsPercStatusQuo_25Quant",
                                                         "CumRealizedGainsPercStatusQuo_50Quant",
                                                         "CumRealizedGainsPercStatusQuo_75Quant",
                                                         "CumRealizedGainsPercStatusQuo_975Quant")

cumGainsTBL <- cumGainsTBL %>%
  bind_cols(cumPotentialGainsPercStatusQuoQuantilesTBL, cumRealizedGainsPercStatusQuoQuantilesTBL)

tail(cumGainsTBL$CumRealizedGainsPercStatusQuo_50Quant, n = 1) #20.3% total realized gains
tail(cumGainsTBL$CumRealizedGainsPercStatusQuo_50Quant[cumGainsTBL$Management %in% c("UnfishedHighBigOld", "UnfishedHighSmallNew")], n = 1) #12.9 of the 20.3% realized gains from high compliance fully protected MPAs
12.9/20.3*100 #i.e., 64% of realized gains is from high compliance fully protected MPAs

tail(cumGainsTBL$CumRealizedGainsPercStatusQuo_025Quant[cumGainsTBL$Management %in% c("UnfishedHighBigOld", "UnfishedHighSmallNew")], n = 1) #14.3
tail(cumGainsTBL$CumRealizedGainsPercStatusQuo_975Quant[cumGainsTBL$Management %in% c("UnfishedHighBigOld", "UnfishedHighSmallNew")], n = 1) #12.0

tail(cumGainsTBL$CumPotentialGainsPercStatusQuo_50Quant, n = 1) #70.8
tail(cumGainsTBL$CumPotentialGainsPercStatusQuo_025Quant, n = 1) #64.8
tail(cumGainsTBL$CumPotentialGainsPercStatusQuo_975Quant, n = 1) #74.7


#Calculate cumulative potential gains for randomly selected sites - to compare with ordered
highCompGainsTBL <- cumGainsTBL %>% 
  dplyr::filter(Management %in% c("UnfishedHighSmallNew", "UnfishedHighBigOld"))

noHighCompGainsTBL <- cumGainsTBL %>% 
  dplyr::filter(!Management %in% c("UnfishedHighSmallNew", "UnfishedHighBigOld"))

randCumPotGainsFilename <- paste0(resultsDir, "RealizedGainsMS_RandomCumulativePotentialGains.rds")
if(file.exists(randCumPotGainsFilename)) {
  randCumPotentialGainsTBL <- readRDS(file = randCumPotGainsFilename)
} else {
  randCumPotentialGainsTBL <- data.frame(matrix(nrow = nrow(cumGainsTBL), ncol = numSamples*numModels)) 
  
  modelNum <- 1
  iterNum <- 0
  for(modelNum in 1:numModels) {
    message("Started model ", modelNum)
    i <- 1
    for (i in 1:numSamples) {
      iterNum <- iterNum + 1
      
      randNoHighRows <- sample(nrow(noHighCompGainsTBL))
      randNoHighTempTBL <- noHighCompGainsTBL[randNoHighRows,]
      randAllTempTBL <- bind_rows(highCompGainsTBL, randNoHighTempTBL)
      
      potentialGainsColname <- paste0("SpaMM_", modelNum, "_PotentialGainsSample_",i)
      statusQuoColname <- paste0("SpaMM_", modelNum, "_StatusQuoPredSample_",i)
      
      randCumPotentialGainsTBL[,iterNum] <- (cumsum(randAllTempTBL[,potentialGainsColname])/sum(randAllTempTBL[,statusQuoColname])*100)[,1]
      message("Random cumulative gains = ", max(randCumPotentialGainsTBL[,iterNum]))
      
    }
    
  }
  #Save this dataset
  saveRDS(object = randCumPotentialGainsTBL, file = randCumPotGainsFilename)
}

randCumPotentialGainsQuantilesTBL <- data.frame(matrix(nrow = nrow(randCumPotentialGainsTBL), ncol = 5))

i = 1
for(i in 1:nrow(randCumPotentialGainsTBL)) {
  randCumPotentialGainsQuantilesTBL[i,] <- unname(quantile(x = randCumPotentialGainsTBL[i,], probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))
}

colnames(randCumPotentialGainsQuantilesTBL) <- c("RandCumPotGainsCompStatusQuo_025Quant",
                                                 "RandCumPotGainsCompStatusQuo_25Quant",
                                                 "RandCumPotGainsCompStatusQuo_50Quant",
                                                 "RandCumPotGainsCompStatusQuo_75Quant",
                                                 "RandCumPotGainsCompStatusQuo_975Quant")

randCumPotentialGainsQuantilesTBL <- randCumPotentialGainsQuantilesTBL %>% 
  bind_cols(cumGainsTBL %>% select(PercentageSites))

#Calculate the median and 95% quantiles for full protection scenario
median(predictTBL$UnfishedHighBigOldBiomassMed) #804 kg/ha
quantile(x = predictTBL$UnfishedHighBigOldBiomassMed, probs = c(0.025, 0.975)) #494 - 1647 kg/ha

##Find out the median and 95% quantiles for expected potential gains at ~30% high compliance
#Get the percentage closest to 30%
closestPerc30 <- randCumPotentialGainsQuantilesTBL$PercentageSites[which.min(abs(randCumPotentialGainsQuantilesTBL$PercentageSites - 30))]

#Maximizing gains
maxPotGains30PercQuant50 <- cumGainsTBL$CumPotentialGainsPercStatusQuo_50Quant[cumGainsTBL$PercentageSites == closestPerc30] #21.8
maxPotGains30PercQuant025 <- cumGainsTBL$CumPotentialGainsPercStatusQuo_025Quant[cumGainsTBL$PercentageSites == closestPerc30] #18.1
maxPotGains30PercQuant975 <- cumGainsTBL$CumPotentialGainsPercStatusQuo_975Quant[cumGainsTBL$PercentageSites == closestPerc30] #24.3

#Random
randPotGains30PercQuant50 <- randCumPotentialGainsQuantilesTBL$RandCumPotGainsCompStatusQuo_50Quant[randCumPotentialGainsQuantilesTBL$PercentageSites == closestPerc30] #14.1
randPotGains30PercQuant025 <- randCumPotentialGainsQuantilesTBL$RandCumPotGainsCompStatusQuo_025Quant[randCumPotentialGainsQuantilesTBL$PercentageSites == closestPerc30] #12.6
randPotGains30PercQuant975 <- randCumPotentialGainsQuantilesTBL$RandCumPotGainsCompStatusQuo_975Quant[randCumPotentialGainsQuantilesTBL$PercentageSites == closestPerc30] #15.7

####>> Fig 3a plot ####
Fig3a_CumPotRealGainsVsPerSites <- ggplot(data = cumGainsTBL,
                                          aes(x = PercentageSites)) +
  #Add a ribbon for the range of values for random siting
  geom_ribbon(data = randCumPotentialGainsQuantilesTBL %>%
                filter(PercentageSites > max(highCompGainsTBL$PercentageSites)),
              aes(ymin = RandCumPotGainsCompStatusQuo_025Quant,
                  ymax = RandCumPotGainsCompStatusQuo_975Quant),
              color = "grey", fill = "grey",
              show.legend = F) +
  #Add a dashed line for the mean random siting
  geom_line(data = randCumPotentialGainsQuantilesTBL %>%
              filter(PercentageSites > max(highCompGainsTBL$PercentageSites)),
            aes(y = RandCumPotGainsCompStatusQuo_50Quant),
            linetype = "dashed", size = 1, color = "black") +
  #Rectangle showing realized and potential gains
  geom_rect(aes(xmin = 100, xmax = 105, ymin = min(CumRealizedGainsPercStatusQuo_025Quant), ymax = 0),
            color = "black", fill = gainColors["RealizedGains"]) +
  geom_rect(aes(xmin = 100, xmax = 105, ymin = 0, ymax = max(CumPotentialGainsPercStatusQuo_975Quant)),
            color = "black", fill = gainColors["PotentialGains"]) +
  geom_segment(x = 100, xend = 105,
               aes(y = max(CumPotentialGainsPercStatusQuo_50Quant),
                   yend = max(CumPotentialGainsPercStatusQuo_50Quant)),
               color = manageColors["UnfishedHighBigOld"], size = 1) +
  geom_segment(x = 100, xend = 105,
               aes(y = min(CumRealizedGainsPercStatusQuo_50Quant),
                   yend = min(CumRealizedGainsPercStatusQuo_50Quant)),
               color = manageColors["Fished"], size = 1) +
  geom_segment(x = 100, xend = 105, y = 0, yend = 0,
               color = "black", size = 2) +
  geom_ribbon(aes(ymin = CumRealizedGainsPercStatusQuo_025Quant,
                  ymax = CumRealizedGainsPercStatusQuo_975Quant),
              fill = manageColors["Fished"],
              alpha = 0.5) +
  geom_line(aes(y = CumRealizedGainsPercStatusQuo_50Quant),
            color = manageColors["Fished"],
            size = 1, linetype = "dashed") +
  geom_ribbon(aes(ymin = CumPotentialGainsPercStatusQuo_025Quant,
                  ymax = CumPotentialGainsPercStatusQuo_975Quant),
              fill = manageColors["UnfishedHighBigOld"],
              alpha = 0.5) +
  geom_line(aes(y = CumPotentialGainsPercStatusQuo_50Quant),
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
                     breaks = c(-20, -10, 0, 10, 20, 30, 40, 50, 60, 70),
                     sec.axis = dup_axis(breaks = c(max(cumGainsTBL$CumPotentialGainsPercStatusQuo_50Quant),
                                                    0,
                                                    min(cumGainsTBL$CumRealizedGainsPercStatusQuo_50Quant)),
                                         labels = c(paste0("Full protection (+", round(max(cumGainsTBL$CumPotentialGainsPercStatusQuo_50Quant)), "%)"),
                                                    "Status quo",
                                                    paste0("All fished (", round(min(cumGainsTBL$CumRealizedGainsPercStatusQuo_50Quant)),"%)")),
                                         name = NULL)) +
  scale_color_manual(values = manageColors,
                     labels = manageLabels) +
  scale_fill_manual(values = manageColors,
                    labels = manageLabels) +
  #Add dotted lines for maximum and random potential gains at 30% protection
  geom_segment(x = closestPerc30, xend = closestPerc30,
               y = 0, yend = maxPotGains30PercQuant50,
               color = "black", linetype = "dotted") +
  geom_segment(x = 0, xend = closestPerc30,
               y = maxPotGains30PercQuant50, yend = maxPotGains30PercQuant50,
               color = "black", linetype = "dotted") +
  geom_segment(x = 0, xend = closestPerc30,
               y = randPotGains30PercQuant50, yend = randPotGains30PercQuant50,
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
            y = max(cumGainsTBL$CumPotentialGainsPercStatusQuo_50Quant)/2,
            label = "Potential",
            color = "white",
            hjust = 0.5,
            vjust = 0.5,
            angle = 270) +
  geom_text(x = 102.5,
            y = min(cumGainsTBL$CumRealizedGainsPercStatusQuo_50Quant)/2,
            label = "Realized",
            color = "black",
            hjust = 0.5,
            vjust = 0.5,
            angle = 270)

######      > b. Map of potential gains ####
###Prepare data
potGainMapTBL <- cumGainsTBL %>%
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

ggsave(plot = Fig3_PotGainsCumBiomassMap, 
       filename = paste0(plotDir, "RealizedGainsMS_Fig3_PotGainsCumGainsMap_", dateNum, ".tiff"),
       width = 11,
       height =  11)

####Calculate the amount of gains from improved compliance alone ####
#Calculate the 95% and 50% quantiles for the cumulative sums from the sample distributions
percGainsLowComplVect <- c()

modelNum <- 1 #for testing
for(modelNum in 1:numModels) {
  message("Started model ", modelNum)
  i = 1
  for(i in 1:numSamples) {
    percGainsLowComplVect <- c(percGainsLowComplVect, sum(predictTBL[predictTBL$Management == "UnfishedLow",paste0("SpaMM_", modelNum, "_PotentialGainsSample_", i)])/sum(predictTBL[,paste0("SpaMM_", modelNum, "_PotentialGainsSample_", i)])*100)
    message("Low compliance gains = ", tail(percGainsLowComplVect, n = 1))
  }
}

lowCompGainsQuants <- quantile(x = percGainsLowComplVect, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))

#Calculate gains from a random selection of fished sites equal to the number of low compliance sites would be
numLowComplSites <- manageSamples$Freq[manageSamples$Var1 == "UnfishedLow"]
propLowComplSites <- numLowComplSites/sum(manageSamples$Freq[manageSamples$Var1 %in% c("UnfishedLow", "UnfishedHighBigOld", "UnfishedHighSmallNew")])
propNumLowComplSites <- round(numLowComplSites*propLowComplSites)
fishedTBL <- predictTBL %>% filter(Management == "Fished")
percGainsFishedRandEqualLowComp <- c()
percGainsFishedRandLowHigh <- c()

modelNum <- 1 #for testing
for(modelNum in 1:numModels) {
  message("Started model ", modelNum)
  i = 1
  for(i in 1:numSamples) {
    #Choose a random selection of fished sites equal to number of low compliance
    randSampleFished <- sample(nrow(fishedTBL), numLowComplSites)
    
    fishedPotGainsIter <- sum(fishedTBL[randSampleFished, paste0("SpaMM_", modelNum, "_PotentialGainsSample_", i)])
    totalPotGainsIter <- sum(predictTBL[,paste0("SpaMM_", modelNum, "_PotentialGainsSample_", i)])
    percGainsFishedRandEqualLowComp <- c(percGainsFishedRandEqualLowComp, fishedPotGainsIter/totalPotGainsIter*100)
    
    #Split the random sample of fished so it matches the average proportion of low compliance sites
    randSampleFishedLow <- sample(randSampleFished, propNumLowComplSites) 
    randSampleFishedHigh <- setdiff(randSampleFished, randSampleFishedLow)
    
    fishedLowPotGainsIter <- sum(fishedTBL[randSampleFishedLow, paste0("SpaMM_", modelNum, "_LowCompGainsSample_", i)])
    fishedHighPotGainsIter <- sum(fishedTBL[randSampleFishedHigh, paste0("SpaMM_", modelNum, "_PotentialGainsSample_", i)])
    
    percGainsFishedRandLowHigh <- c(percGainsFishedRandLowHigh, (fishedLowPotGainsIter + fishedHighPotGainsIter)/totalPotGainsIter*100)
    message("Fished gains = ", tail(percGainsFishedRandEqualLowComp, n = 1), " (perfect compl); ",
            tail(percGainsFishedRandLowHigh, n = 1), " (average compl)")
  }
}

fishedGainsPerfComplQuants <- quantile(x = percGainsFishedRandEqualLowComp, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
fishedGainsAvgComplQuants <- quantile(x = percGainsFishedRandLowHigh, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))