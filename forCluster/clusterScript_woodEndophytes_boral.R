# script to fit boral models on the cluster
# saves each model run as a separate file

# file structure:
# code/boral_modfit.R
# code/load_data.R
# derived_data/complete_subset_list.RData
# output

# load libraries -----------------------------------------------------
library(dplyr)
library(tidyr)
library(boral)

# load functions -----------------------------------------------------
source("code/boral_modfit.R")
source("code/load_data.R")

# load data -----------------------------------------------------
complete_subset_list <- readRDS(file = "derived_data/complete_subset_list.RData")
runType <- "real"


# LV-only model -----------------------------------------------------
NUMRUNS <- 12
SEED <- .Random.seed[1:NUMRUNS] #set seed for different runs
for(i in 1:NUMRUNS){
    curr.fit.stats <- fit_boral_lvOnly(complete_subset_list, seed = SEED[i], runType = runType)
    curr.fileName <- paste0("output/LV-only_run", i, "_runType", runType, ".RData")
    saveRDS(curr.fit.stats, file = curr.fileName)
}

# Traits-and-LVs model (all X vars) -----------------------------------------------------
NUMRUNS <- 12
SEED <- .Random.seed[1:NUMRUNS] #set seed for different runs
for(i in 1:NUMRUNS){
    curr.fit.stats <- fit_boral_allX(complete_subset_list, seed = SEED[i], runType = runType)
    curr.fileName <- paste0("output/Traits-and-LVs_allX_run", i, "_runType", runType, ".RData")
    saveRDS(curr.fit.stats, file = curr.fileName)
}

# Traits-and-LVs model (select X vars) -----------------------------------------------------
NUMRUNS <- 12
SEED <- .Random.seed[1:NUMRUNS] #set seed for different runs
for(i in 1:NUMRUNS){
    curr.fit.stats <- fit_boral_selectX(complete_subset_list, seed = SEED[i], runType = runType)
    curr.fileName <- paste0("output/Traits-and-LVs_selectX_run", i, "_runType", runType, ".RData")
    saveRDS(curr.fit.stats, file = curr.fileName)
    }
