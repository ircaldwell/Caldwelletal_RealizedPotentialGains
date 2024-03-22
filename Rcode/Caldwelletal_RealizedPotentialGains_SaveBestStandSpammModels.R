###### Caldwell et al. Code for "Protection efforts have resulted in ~10% of existing fish biomass on global coral reefs" ####
######  Runs and saves the best spaMM models  ####
######      Author: Iain R. Caldwell
######      Last revised: Mar 21, 2024
######    1. Open the files with the best spaMM results (and filenames), and the fish survey data with covariates for running the models
######    2. Change the formulas so all continuous predictors are standardized
######    3. Run and save the standardized best models
rm(list = ls()) #remove past stored objects

####Load packages and libraries ####
library(tidyverse)
library(spaMM) #to run spatial models

####  Set parameters and directories ####
dataDir <- "../Data/"
modelDir <- "../Models/"

######    1. Open the files with the best spaMM results (and filenames), then process the fish survey data for re-running models ####
bestSpammResTBL <- read_csv(paste0(dataDir, "Caldwelletal_RealizedPotentialGains_BestSpammSummResWithFilenames.csv")) 

# Open file wih fish surveys and covariates
fishSurveyTBL <- read_csv(file = paste0(dataDir, "Caldwelletal_RealizedPotentialGains_FishSurveysCovariates.csv"))

######    2. Change the formulas so all continuous predictors are standardized ####
#### Set up continuous columns ####
gravCols <- c("logGrav_NearMarket", "logGrav_NearPop")
sstCols <- c("SST_mean_2yr", "SST_max_2yr", "SST_sd_2yr", "logSST_kurtosis_2yr", "SST_skewness_2yr")
sstaCols <- c("SSTa_mean_2yr") 
dhwCols <- c("logDHW_mean_2yr", "logDHW_max_2yr")
parCols <- c("PAR_mean_2yr", "PAR_sd_2yr", "PAR_skewness_2yr", "logPAR_kurtosis_2yr")
chlaCols <- c("logChlA_mean_2yr", "logChlA_max_2yr")
msecCols <- c("logwave_mean")
continuousPreds <- c(gravCols, sstCols, sstaCols, dhwCols, chlaCols, parCols, msecCols)

contPred <- continuousPreds[1]
for(contPred in continuousPreds) {
  message("Added scale to ", contPred)
  bestSpammResTBL <- bestSpammResTBL %>% 
    mutate(Formula = gsub(pattern = contPred, replacement = paste0("scale(", contPred, ")"), x = Formula, fixed = T))
}

######    3. Run and save the standardized best models ####
# Function to run and save the models
saveSpammModelsFun <- function(i, spammFormulaDataset, spammDataset, saveDir) {
  #Run spaMM model
  spammModel <- spaMM::fitme(formula = as.formula(spammFormulaDataset$Formula[i]),
                             data = spammDataset,
                             family = Gamma(log))
  
  #Save the spaMM model
  saveRDS(object = spammModel, file = paste0(saveDir, spammFormulaDataset$ModelFilename[i]))
}

#Check if the model folder exists and create it if not
if(!file.exists(modelDir)) {
  dir.create(file.path(modelDir))
} 

#Check if the models have already been run and saved
existingModelNames <- list.files(path = modelDir)
remModelsToRunTBL <- bestSpammResTBL %>% 
  filter(!ModelFilename %in% existingModelNames)

#Run any remaining models
if(nrow(remModelsToRunTBL) > 0) {
  iter = 1:nrow(remModelsToRunTBL)
  n.cores = min(c(length(iter), detectCores()-1))
  
  start_time <- Sys.time()
  message("Saving best spaMM models - started at ", start_time) 
  
  clust = makeCluster(n.cores)
  
  clusterExport(cl = clust, varlist = c("fishSurveyTBL", "remModelsToRunTBL", "saveSpammModelsFun", "fitme", "tibble"))
  parSapply(cl = clust, X = iter, FUN = saveSpammModelsFun,
            spammFormulaDataset = remModelsToRunTBL, spammDataset = fishSurveyTBL, saveDir = modelDir)
  
  stopCluster(cl = clust)
  
  end_time <- Sys.time()
  message("Finished saving best spaMM models - ended at ", end_time) 
  
  elapsedTime = end_time - start_time
  message("spaMM models took ", elapsedTime, " hours")
} else {
  message("All the best spaMM models have been run")
}

