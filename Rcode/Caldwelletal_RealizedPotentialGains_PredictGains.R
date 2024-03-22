###### Caldwell et al. Code for "Protection efforts have resulted in ~10% of existing fish biomass on global coral reefs" ####
######  Predicts coral reef fish biomass under counterfactual scenarios and calculates fish biomass gains using best spaMM models  ####
######      Code author: Iain R. Caldwell
######      Last revised: March 21, 2024
######    1. Open the file with the best biomass spaMM results (Gamma) and the survey data
######    2. Loop through all of the best models and get predictions with error
######        a) Status quo biomass (standardizing methods, etc) 
######        b) All fished scenario biomass
######        c) All fully protected scenario biomass
######        d) All restricted fishing biomass
######    3. Calculate multi-model means with combined between and within model error
######        a) Status quo biomass (standardizing methods, etc) 
######        b) All fished scenario biomass
######        c) All fully protected scenario biomass
######        d) All restricted fishing biomass
######    4. Simulate distributions for each of the scenarios and use those to calculate back transformed gains with uncertainties    
######        a) Realized gains (Status Quo - Fished)  
######        b) Potential gains (Unfished High Big Old scenario - Status Quo) 
######        c) Restricted fishing gains (Restricted scenario - Status Quo) 
######    5. Use simulated distributions to calculate biomasses and gains as a percentage of status quo biomass     
######        a) Status quo biomass as a percentage of total status quo biomass (For Fig 1a)
######        b) Fished biomass as a percentage of total status quo biomass (for Fig 1b)
######        c) Realized gains (Status Quo - Fished) % of total status quo 
######        d) Potential gains (Unfished High Big Old scenario - Status Quo) as % of status quo
######        e) Restricted fishing gains (Restricted scenario - Status Quo) as % of status quo
######    6. Calculate median predictions with 95% quantiles from all sampled distributions ####
######        a) Status quo (biomass and % of total status quo)  
######        b) all fished scenario (biomass and % of total status quo)
######        c) all fully protected scenario (biomass)
######        d) restricted scenario biomass (biomass)
######        e) realized gains (biomass and % of total status quo)
######        f) potential gains (biomass and % of total status quo)
######        g) restricted gains (biomass and % of total status quo)
######    7. Save results file
rm(list = ls()) #remove past stored objects
options(scipen = 999) #turn off scientific notation
set.seed(1234)

####Load packages and libraries ####
library(tidyverse)
library(spaMM) #to run spatial models
library(plotrix) #to calculate SE among models

####  Set parameters and directories ####
numSamples <- 1000
dataDir <- "../Data/"
modelDir <- "../Models/"

######    1. Open the file with the best biomass spaMM results (Gamma) and the survey data ####
bestSpammResTBL <- read_csv(paste0(dataDir, "Caldwelletal_RealizedPotentialGains_BestSpammSummResWithFilenames.csv")) 

#Check if the next 2 sections have already been run and only run if not (since it takes a long time)
fishPredMeanSDfilename <- paste0(dataDir, "Caldwelletal_RealizedPotentialGains_FishBiomassPredMeanSD.rds")

if(file.exists(fishPredMeanSDfilename)) {
  fishSurveyTBL <- readRDS(fishPredMeanSDfilename)
} else {
  
  # Open file wih fish surveys and covariates
  fishSurveyTBL <- read_csv(file = paste0(dataDir, "Caldwelletal_RealizedPotentialGains_FishSurveysCovariates.csv"))
  
  ######    2. Loop through all of the best models and get predictions with error ####
  ######        a) Status quo biomass (standardizing methods, etc) 
  ######        b) All fished scenario biomass
  ######        c) All fully protected scenario biomass
  ######        d) All restricted fishing biomass
  
  statusQuoTBL <- fishSurveyTBL %>% 
    mutate(Habitat = "Slope",
           Depth = "4-10m",
           CensusMethod = "Belt transect")
  
  allFishedTBL <- statusQuoTBL %>% 
    mutate(Management = "Fished")
  
  allProtTBL <- statusQuoTBL %>% 
    mutate(Management = "UnfishedHighBigOld")
  
  allRestTBL <- statusQuoTBL %>% 
    mutate(Management = "Restricted")
  
  ### Predicting values in log space so that they are same scale as SD, then convert back to biomass to calculate gains
  i = 1
  for(i in 1:nrow(bestSpammResTBL)) {
    message("Started model ", i, " of ", nrow(bestSpammResTBL))
    iterSpaMM <- readRDS(file = paste0(modelDir, bestSpammResTBL$ModelFilename[i]))
    
    modelPredTBL <- tibble(StatusQuoLogBiomassPred = unname(predict(object = iterSpaMM,
                                                                    type = "link", #Predict on log scale for later calculating averages
                                                                    binding = NA,
                                                                    newdata = statusQuoTBL,
                                                                    re.form = as.formula(gsub(pattern = "+Matern(1|Easting+Northing)",
                                                                                              replacement = "",
                                                                                              x = bestSpammResTBL$Formula[i],
                                                                                              fixed = T)))),
                           #Calculating SE from point prediction variance
                           StatusQuoLogBiomassSD = sqrt(unname(get_predVar(object = iterSpaMM,
                                                                           newdata = statusQuoTBL, 
                                                                           re.form = as.formula(gsub(pattern = "+Matern(1|Easting+Northing)",
                                                                                                     replacement = "",
                                                                                                     x = bestSpammResTBL$Formula[i],
                                                                                                     fixed = T))))),
                           FishedLogBiomassPred = unname(predict(object = iterSpaMM,
                                                                 type = "link", #Predict on log scale for later calculating averages
                                                                 binding = NA,
                                                                 newdata = allFishedTBL,
                                                                 re.form = as.formula(gsub(pattern = "+Matern(1|Easting+Northing)",
                                                                                           replacement = "",
                                                                                           x = bestSpammResTBL$Formula[i],
                                                                                           fixed = T)))),
                           FishedLogBiomassSD = sqrt(unname(get_predVar(object = iterSpaMM,
                                                                        newdata = allFishedTBL, 
                                                                        re.form = as.formula(gsub(pattern = "+Matern(1|Easting+Northing)",
                                                                                                  replacement = "",
                                                                                                  x = bestSpammResTBL$Formula[i],
                                                                                                  fixed = T))))),
                           UnfishedHighBigOldLogBiomassPred = unname(predict(object = iterSpaMM,
                                                                             type = "link", #Predict on log scale for later calculating averages
                                                                             binding = NA,
                                                                             newdata = allProtTBL,
                                                                             re.form = as.formula(gsub(pattern = "+Matern(1|Easting+Northing)",
                                                                                                       replacement = "",
                                                                                                       x = bestSpammResTBL$Formula[i],
                                                                                                       fixed = T)))),
                           UnfishedHighBigOldLogBiomassSD = sqrt(unname(get_predVar(object = iterSpaMM,
                                                                                    newdata = allProtTBL, 
                                                                                    re.form = as.formula(gsub(pattern = "+Matern(1|Easting+Northing)",
                                                                                                              replacement = "",
                                                                                                              x = bestSpammResTBL$Formula[i],
                                                                                                              fixed = T))))),
                           RestrictedLogBiomassPred = unname(predict(object = iterSpaMM,
                                                                     type = "link", #Predict on log scale for later calculating averages
                                                                     binding = NA,
                                                                     newdata = allRestTBL,
                                                                     re.form = as.formula(gsub(pattern = "+Matern(1|Easting+Northing)",
                                                                                               replacement = "",
                                                                                               x = bestSpammResTBL$Formula[i],
                                                                                               fixed = T)))),
                           RestrictedLogBiomassSD = sqrt(unname(get_predVar(object = iterSpaMM,
                                                                            newdata = allRestTBL, 
                                                                            re.form = as.formula(gsub(pattern = "+Matern(1|Easting+Northing)",
                                                                                                      replacement = "",
                                                                                                      x = bestSpammResTBL$Formula[i],
                                                                                                      fixed = T)))))) 
    
    #Truncate realized and potential gains at zero
    modelPredTBL <- modelPredTBL %>% 
      mutate(FishedLogBiomassPred = ifelse(FishedLogBiomassPred > StatusQuoLogBiomassPred, StatusQuoLogBiomassPred, FishedLogBiomassPred), #truncate at 0
             UnfishedHighBigOldLogBiomassPred = ifelse(StatusQuoLogBiomassPred > UnfishedHighBigOldLogBiomassPred, StatusQuoLogBiomassPred, UnfishedHighBigOldLogBiomassPred))
    
    #Rename columns with the model number
    colnames(modelPredTBL) <- paste0("SpaMM_", i, "_", colnames(modelPredTBL))
    
    fishSurveyTBL <- fishSurveyTBL %>% 
      bind_cols(modelPredTBL)
    
  }
  
  ######    3. Calculate multi-model means with combined between and within model error ####
  ######        a) Status quo biomass (standardizing methods, etc) 
  ######        b) All fished scenario biomass
  ######        c) All fully protected scenario biomass
  ######        d) All restricted fishing biomass
  
  fishSurveyTBL <- fishSurveyTBL %>%
    mutate(StatusQuoLogBiomassMean = rowMeans(select(., contains("_StatusQuoLogBiomassPred"))),
           StatusQuoLogBiomassWithinModelSD = sqrt(rowMeans(select(., contains("_StatusQuoLogBiomassSD"))^2)),
           FishedLogBiomassMean = rowMeans(select(., contains("_FishedLogBiomassPred"))),
           FishedLogBiomassWithinModelSD = sqrt(rowMeans(select(., contains("_FishedLogBiomassSD"))^2)),
           UnfishedHighBigOldLogBiomassMean = rowMeans(select(., contains("_UnfishedHighBigOldLogBiomassPred"))),
           UnfishedHighBigOldLogBiomassWithinModelSD = sqrt(rowMeans(select(., contains("_UnfishedHighBigOldLogBiomassSD"))^2)),
           RestrictedLogBiomassMean = rowMeans(select(., contains("_RestrictedLogBiomassPred"))),
           RestrictedLogBiomassWithinModelSD = rowMeans(select(., contains("_RestrictedLogBiomassSD"))^2))
  
  #Also calculate the standard error around the means (from variation among models) and combine with the propagated uncertainty
  fishSurveyTBL <- fishSurveyTBL %>% 
    rowwise() %>% 
    mutate(StatusQuoLogBiomassAmongModelSD = std.error(c_across(contains("_StatusQuoLogBiomassPred")))*sqrt(nrow(bestSpammResTBL)),
           FishedLogBiomassAmongModelSD = std.error(c_across(contains("_FishedLogBiomassPred")))*sqrt(nrow(bestSpammResTBL)),
           UnfishedHighBigOldLogBiomassAmongModelSD = std.error(c_across(contains("_UnfishedHighBigOldLogBiomassPred")))*sqrt(nrow(bestSpammResTBL)),
           RestrictedLogBiomassAmongModelSD = std.error(c_across(contains("_RestrictedLogBiomassPred")))*sqrt(nrow(bestSpammResTBL))) %>% 
    ungroup() %>% 
    #Calculate the total uncertainty as the sum in quadrature of the within and among model error
    mutate(StatusQuoLogBiomassSD = sqrt((StatusQuoLogBiomassWithinModelSD^2) + (StatusQuoLogBiomassAmongModelSD^2)),
           FishedLogBiomassSD = sqrt((FishedLogBiomassWithinModelSD^2) + (FishedLogBiomassAmongModelSD^2)),
           UnfishedHighBigOldLogBiomassSD = sqrt((UnfishedHighBigOldLogBiomassWithinModelSD^2) + (UnfishedHighBigOldLogBiomassAmongModelSD^2)),
           RestrictedLogBiomassSD = sqrt((RestrictedLogBiomassWithinModelSD^2) + (RestrictedLogBiomassAmongModelSD^2)))
  
  #Save this file so it can be called instead of re-running this (since it takes a long time)
  saveRDS(object = fishSurveyTBL,
          file = fishPredMeanSDfilename)
  
}

######    4. Simulate distributions for each of the scenarios and use those to calculate back transformed gains with uncertainties ####     
######        a) Realized gains (Status Quo - Fished)  
######        b) Potential gains (Unfished High Big Old scenario - Status Quo) 
######        c) Restricted fishing gains (Restricted scenario - Status Quo) 

#Check if the next section has already been run and only run if not (since it also takes a long time)
fishPredWithSimsfilename <- paste0(dataDir, "Caldwelletal_RealizedPotentialGains_FishBiomassPredSims.rds")

if(file.exists(fishPredWithSimsfilename)) {
  fishSurveyTBL <- readRDS(file = fishPredWithSimsfilename)
} else {
  #Create blank dataframes to fill in samples from each of the distributions
  statusQuoSampleDist <- data.frame(matrix(ncol = numSamples, nrow = nrow(fishSurveyTBL)))
  colnames(statusQuoSampleDist) <- gsub(pattern = "X",
                                        replacement = "StatusQuoPredSample_",
                                        x = colnames(statusQuoSampleDist))
  
  fishedSampleDist <- data.frame(matrix(ncol = numSamples, nrow = nrow(fishSurveyTBL)))
  colnames(fishedSampleDist) <- gsub(pattern = "X",
                                     replacement = "FishedPredSample_",
                                     x = colnames(fishedSampleDist))
  
  unfishedSampleDist <- data.frame(matrix(ncol = numSamples, nrow = nrow(fishSurveyTBL)))
  colnames(unfishedSampleDist) <- gsub(pattern = "X",
                                       replacement = "UnfishedPredSample_",
                                       x = colnames(unfishedSampleDist))
  
  restSampleDist <- data.frame(matrix(ncol = numSamples, nrow = nrow(fishSurveyTBL)))
  colnames(restSampleDist) <- gsub(pattern = "X",
                                   replacement = "RestPredSample_",
                                   x = colnames(restSampleDist))
  
  realizedGainsSampleDist <- data.frame(matrix(ncol = numSamples, nrow = nrow(fishSurveyTBL))) 
  colnames(realizedGainsSampleDist) <- gsub(pattern = "X",
                                            replacement = "RealizedGainsSample_",
                                            x = colnames(realizedGainsSampleDist))
  
  potentialGainsSampleDist <- data.frame(matrix(ncol = numSamples, nrow = nrow(fishSurveyTBL))) 
  colnames(potentialGainsSampleDist) <- gsub(pattern = "X",
                                             replacement = "PotentialGainsSample_",
                                             x = colnames(potentialGainsSampleDist))
  
  
  restGainsSampleDist <- data.frame(matrix(ncol = numSamples, nrow = nrow(fishSurveyTBL)))
  colnames(restGainsSampleDist) <- gsub(pattern = "X",
                                        replacement = "RestGainsSample_",
                                        x = colnames(restGainsSampleDist))
  
  #Fill in dataframes with samples frpm distributions
  j = 1
  for(j in 1:nrow(fishSurveyTBL)) {
    message("Started row ", j, " of ", nrow(fishSurveyTBL))
    iterManagement <- fishSurveyTBL$Management[j]
    
    statusQuoSampleDist[j,] <- unname(exp(rnorm(n = numSamples,
                                                mean = fishSurveyTBL$StatusQuoLogBiomassMean[j],
                                                sd = fishSurveyTBL$StatusQuoLogBiomassSD[j])))
    
    if(iterManagement == "Fished") {
      fishedSampleDist[j,] <- statusQuoSampleDist[j,]
    } else {
      fishedSampleDist[j,] <- unname(exp(rnorm(n = numSamples,
                                               mean = fishSurveyTBL$FishedLogBiomassMean[j],
                                               sd = fishSurveyTBL$FishedLogBiomassSD[j])))
      fishedSampleDist[j, fishedSampleDist[j,] > statusQuoSampleDist[j,]] <-
        statusQuoSampleDist[j, fishedSampleDist[j,] > statusQuoSampleDist[j,]]
    }
    
    if(iterManagement == "UnfishedHighBigOld") {
      unfishedSampleDist[j,] <- statusQuoSampleDist[j,]
    } else {
      unfishedSampleDist[j,] <- unname(exp(rnorm(n = numSamples,
                                                 mean = fishSurveyTBL$UnfishedHighBigOldLogBiomassMean[j],
                                                 sd = fishSurveyTBL$UnfishedHighBigOldLogBiomassSD[j])))
      unfishedSampleDist[j, unfishedSampleDist[j,] < statusQuoSampleDist[j,]] <-
        statusQuoSampleDist[j, unfishedSampleDist[j,] < statusQuoSampleDist[j,]]
    }
    
    if(iterManagement == "Restricted") {
      restSampleDist[j,] <- statusQuoSampleDist[j,]
    } else {
      restSampleDist[j,] <- unname(exp(rnorm(n = numSamples,
                                             mean = fishSurveyTBL$RestrictedLogBiomassMean[j],
                                             sd = fishSurveyTBL$RestrictedLogBiomassSD[j])))
    }
    
    if(iterManagement == "Fished") {
      realizedGainsSampleDist[j,] <- 0
    } else {
      realizedGainsSampleDist[j,] <- statusQuoSampleDist[j,] - fishedSampleDist[j,]
      #realizedGainsSampleDist[j, realizedGainsSampleDist[j,] < 0] <- 0 #Remove any negative realized gains
    }
    
    if(iterManagement == "UnfishedHighBigOld") {
      potentialGainsSampleDist[j,] <- 0
    } else {
      potentialGainsSampleDist[j,] <- unfishedSampleDist[j,] - statusQuoSampleDist[j,]
      #potentialGainsSampleDist[j, potentialGainsSampleDist[j,] < 0] <- 0 #Remove any negative potential gains
    }
    
    if(iterManagement == "Restricted") {
      restGainsSampleDist[j,] <- 0
    } else {
      restGainsSampleDist[j,] <- restSampleDist[j,] - statusQuoSampleDist[j,]
    }
    
  }
  
  allModelSimsTBL <- bind_cols(statusQuoSampleDist, fishedSampleDist, unfishedSampleDist, restSampleDist,
                               realizedGainsSampleDist, potentialGainsSampleDist, restGainsSampleDist)
  
  fishSurveyTBL <- fishSurveyTBL %>% 
    bind_cols(allModelSimsTBL)
  
  saveRDS(object = fishSurveyTBL, fishPredWithSimsfilename)
  
  #Check that all fished samples are less than or equal to status quo and all unfished samples are greater than or equal to status quo and there are no zeroes in the realized or potential gains
  sum(fishSurveyTBL$FishedPredSample_1 > fishSurveyTBL$StatusQuoPredSample_1) 
  range(fishSurveyTBL %>% select(contains("RealizedGainsSample"))) 
  sum(fishSurveyTBL$UnfishedPredSample_1 < fishSurveyTBL$StatusQuoPredSample_1)
  range(fishSurveyTBL %>% select(contains("PotentialGainsSample")))  
}

######    5. Use simulated distributions to calculate biomasses and gains as a percentage of status quo biomass ####     
######        a) Status quo biomass as a percentage of total status quo biomass (For Fig 1a)
######        b) Fished biomass as a percentage of total status quo biomass (for Fig 1b)
######        c) Realized gains (Status Quo - Fished) % of total status quo 
######        d) Potential gains (Unfished High Big Old scenario - Status Quo) as % of status quo
######        e) Restricted fishing gains (Restricted scenario - Status Quo) as % of status quo

#Check if the next section has already been run and only run if not (since it also takes a long time)
fishPredAllSimsfilename <- paste0(dataDir, "Caldwelletal_RealizedPotentialGains_FishBiomassPredAllSims.rds")

if(file.exists(fishPredAllSimsfilename)) {
  fishSurveyTBL <- readRDS(fishPredAllSimsfilename)
} else {
  #Create new empty dataframes to fill
  statusQuoPercTotalStatusQuoSampleDist <- data.frame(matrix(ncol = numSamples, nrow = nrow(fishSurveyTBL)))
  colnames(statusQuoPercTotalStatusQuoSampleDist) <- gsub(pattern = "X",
                                                          replacement = "StatusQuoPercTotalStatusQuoSample_",
                                                          x = colnames(statusQuoPercTotalStatusQuoSampleDist))
  
  fishedPercTotalStatusQuoSampleDist <- data.frame(matrix(ncol = numSamples, nrow = nrow(fishSurveyTBL)))
  colnames(fishedPercTotalStatusQuoSampleDist) <- gsub(pattern = "X",
                                                       replacement = "FishedPercTotalStatusQuoSample_",
                                                       x = colnames(fishedPercTotalStatusQuoSampleDist))
  
  realizedGainsPercStatusQuoSampleDist <- data.frame(matrix(ncol = numSamples, nrow = nrow(fishSurveyTBL))) 
  colnames(realizedGainsPercStatusQuoSampleDist) <- gsub(pattern = "X",
                                                         replacement = "RealizedGainsPercStatusQuoSample_",
                                                         x = colnames(realizedGainsPercStatusQuoSampleDist))
  
  potentialGainsPercStatusQuoSampleDist <- data.frame(matrix(ncol = numSamples, nrow = nrow(fishSurveyTBL))) 
  colnames(potentialGainsPercStatusQuoSampleDist) <- gsub(pattern = "X",
                                                          replacement = "PotentialGainsPercStatusQuoSample_",
                                                          x = colnames(potentialGainsPercStatusQuoSampleDist))
  
  restGainsPercStatusQuoSampleDist <- data.frame(matrix(ncol = numSamples, nrow = nrow(fishSurveyTBL)))
  colnames(restGainsPercStatusQuoSampleDist) <- gsub(pattern = "X",
                                                     replacement = "RestGainsPercStatusQuoSample_",
                                                     x = colnames(restGainsPercStatusQuoSampleDist))
  #Calculate the total status quo biomass
  totalStatusQuoBiomassSampleDist <- fishSurveyTBL %>% 
    summarise(across(contains("StatusQuoPredSample_"), sum)) 
  
  statusQuoSampleDist <- fishSurveyTBL %>% 
    select(starts_with("StatusQuoPredSample_"))
  
  fishedSampleDist <- fishSurveyTBL %>% 
    select(starts_with("FishedPredSample_"))
  
  realizedGainsSampleDist <- fishSurveyTBL %>% 
    select(starts_with("RealizedGainsSample_"))
  
  potentialGainsSampleDist <- fishSurveyTBL %>% 
    select(starts_with("PotentialGainsSample_"))
  
  restGainsSampleDist <- fishSurveyTBL %>% 
    select(starts_with("RestGainsSample_"))
  
  i = 1 #for testing
  for(i in 1:nrow(fishSurveyTBL)) {
    message("Started row ", i, " of ", nrow(fishSurveyTBL))
    iterManagement <- fishSurveyTBL$Management[i]
    
    statusQuoPercTotalStatusQuoSampleDist[i,] <- statusQuoSampleDist[i,]/totalStatusQuoBiomassSampleDist*100
    fishedPercTotalStatusQuoSampleDist[i,] <- fishedSampleDist[i,]/totalStatusQuoBiomassSampleDist*100
    realizedGainsPercStatusQuoSampleDist[i,] <- realizedGainsSampleDist[i,]/totalStatusQuoBiomassSampleDist*100
    potentialGainsPercStatusQuoSampleDist[i,] <- potentialGainsSampleDist[i,]/totalStatusQuoBiomassSampleDist*100
    restGainsPercStatusQuoSampleDist[i,] <- restGainsSampleDist[i,]/totalStatusQuoBiomassSampleDist*100
  }
  
  #Check that all ranges make sense
  range(statusQuoPercTotalStatusQuoSampleDist) 
  range(fishedPercTotalStatusQuoSampleDist) 
  range(realizedGainsPercStatusQuoSampleDist)
  range(potentialGainsPercStatusQuoSampleDist)
  range(restGainsPercStatusQuoSampleDist)
  
  #Combine these files with the survey and covariate data and save again
  allGainSimsTBL <- bind_cols(statusQuoPercTotalStatusQuoSampleDist,
                              fishedPercTotalStatusQuoSampleDist,
                              realizedGainsPercStatusQuoSampleDist,
                              potentialGainsPercStatusQuoSampleDist,
                              restGainsPercStatusQuoSampleDist)
  
  fishSurveyTBL <- fishSurveyTBL %>% 
    bind_cols(allGainSimsTBL)
  
  saveRDS(object = fishSurveyTBL, fishPredAllSimsfilename)
}

######    6. Calculate median predictions with 95% quantiles from all sampled distributions ####
######        a) Status quo (biomass and % of total status quo)  
######        b) all fished scenario (biomass and % of total status quo)
######        c) all fully protected scenario (biomass)
######        d) restricted scenario biomass (biomass)
######        e) realized gains (biomass and % of total status quo)
######        f) potential gains (biomass and % of total status quo)
######        g) restricted gains (biomass and % of total status quo)

fishSurveyTBL <- fishSurveyTBL %>% 
  rowwise() %>% 
  mutate(StatusQuoBiomassMed = median(c_across(starts_with("StatusQuoPredSample_"))),
         StatusQuoBiomassLow025 = unname(quantile(c_across(starts_with("StatusQuoPredSample_")), probs = 0.025)),
         StatusQuoBiomassHigh975 = unname(quantile(c_across(starts_with("StatusQuoPredSample_")), probs = 0.975)),
         StatusQuoPercTotalStatusQuoMed = median(c_across(starts_with("StatusQuoPercTotalStatusQuoSample_"))),
         StatusQuoPercTotalStatusQuoLow025 = unname(quantile(c_across(starts_with("StatusQuoPercTotalStatusQuoSample_")), probs = 0.025)),
         StatusQuoPercTotalStatusQuoLow975 = unname(quantile(c_across(starts_with("StatusQuoPercTotalStatusQuoSample_")), probs = 0.975)),
         FishedBiomassMed = median(c_across(starts_with("FishedPredSample_"))),
         FishedBiomassLow025 = unname(quantile(c_across(starts_with("FishedPredSample_")), probs = 0.025)),
         FishedBiomassHigh975 = unname(quantile(c_across(starts_with("FishedPredSample_")), probs = 0.975)),
         FishedPercTotalStatusQuoMed = median(c_across(starts_with("FishedPercTotalStatusQuoSample_"))),
         FishedPercTotalStatusQuoLow025 = unname(quantile(c_across(starts_with("FishedPercTotalStatusQuoSample_")), probs = 0.025)),
         FishedPercTotalStatusQuoLow975 = unname(quantile(c_across(starts_with("FishedPercTotalStatusQuoSample_")), probs = 0.975)),
         UnfishedHighBigOldBiomassMed = median(c_across(starts_with("UnfishedPredSample_"))),
         UnfishedHighBigOldBiomassLow025 = unname(quantile(c_across(starts_with("UnfishedPredSample_")), probs = 0.025)),
         UnfishedHighBigOldBiomassHigh975 = unname(quantile(c_across(starts_with("UnfishedPredSample_")), probs = 0.975)),
         RestBiomassMed = median(c_across(starts_with("RestPredSample_"))),
         RestBiomassLow025 = unname(quantile(c_across(starts_with("RestPredSample_")), probs = 0.025)),
         RestBiomassHigh975 = unname(quantile(c_across(starts_with("RestPredSample_")), probs = 0.975)),
         RealizedGainMed = median(c_across(starts_with("RealizedGainsSample_"))),
         RealizedGainLow025 = unname(quantile(c_across(starts_with("RealizedGainsSample_")), probs = 0.025)),
         RealizedGainHigh975 = unname(quantile(c_across(starts_with("RealizedGainsSample_")), probs = 0.975)),
         RealizedGainPercTotalStatusQuoMed = median(c_across(starts_with("RealizedGainsPercStatusQuoSample_"))),
         RealizedGainPercTotalStatusQuoLow025 = unname(quantile(c_across(starts_with("RealizedGainsPercStatusQuoSample_")), probs = 0.025)),
         RealizedGainPercTotalStatusQuoHigh975 = unname(quantile(c_across(starts_with("RealizedGainsPercStatusQuoSample_")), probs = 0.975)),
         PotentialGainMed = median(c_across(starts_with("PotentialGainsSample_"))),
         PotentialGainLow025 = unname(quantile(c_across(starts_with("PotentialGainsSample_")), probs = 0.025)),
         PotentialGainHigh975 = unname(quantile(c_across(starts_with("PotentialGainsSample_")), probs = 0.975)),
         PotentialGainPercTotalStatusQuoMed = median(c_across(starts_with("PotentialGainsPercStatusQuoSample_"))),
         PotentialGainPercTotalStatusQuoLow025 = unname(quantile(c_across(starts_with("PotentialGainsPercStatusQuoSample_")), probs = 0.025)),
         PotentialGainPercTotalStatusQuoHigh975 = unname(quantile(c_across(starts_with("PotentialGainsPercStatusQuoSample_")), probs = 0.975)),
         RestGainMed = median(c_across(starts_with("RestGainsSample_"))),
         RestGainLow025 = unname(quantile(c_across(starts_with("RestGainsSample_")), probs = 0.025)),
         RestGainHigh975 = unname(quantile(c_across(starts_with("RestGainsSample_")), probs = 0.975)),
         RestGainPercTotalStatusQuoMed = median(c_across(starts_with("RestGainsPercStatusQuoSample_"))),
         RestGainPercTotalStatusQuoLow025 = unname(quantile(c_across(starts_with("RestGainsPercStatusQuoSample_")), probs = 0.025)),
         RestGainPercTotalStatusQuoHigh975 = unname(quantile(c_across(starts_with("RestGainsPercStatusQuoSample_")), probs = 0.975))) %>% 
  ungroup()

######    7. Save results file ####
spammPredFilename = paste0(dataDir, "Caldwelletal_RealizedPotentialGains_FishBiomassGainsPredSimsCovariates.rds")
saveRDS(object = fishSurveyTBL, file = spammPredFilename)

#Rename to predictTBL and remove the fishSurveyTBL in case it is run through source
predictTBL <- fishSurveyTBL
rm(fishSurveyTBL)
