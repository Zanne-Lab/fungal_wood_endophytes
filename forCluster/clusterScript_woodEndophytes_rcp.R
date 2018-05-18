# script to fit rcp models on the cluster

# file structure:
# code/rcp_modfit.R
# derived_data/complete_subset_list.RData
# output

# load libraries -----------------------------------------------------
library(dplyr)
library(tidyr)
library(methods) # workaround so functions can be found when running via RScript
library(Matrix) # sameo
library(glmnet) # sameo
library(RCPmod) # normally, this is the only package you'll need to load, rest will load via namespace


# load functions -----------------------------------------------------
source("code/rcp_modfit.R")

# load data -----------------------------------------------------
complete_subset_list <- readRDS(file = "derived_data/complete_subset_list.RData")


# specify a vector of nRCP values -----------------------------------------------------
nRCP.lower <- 2
nRCP.upper <- 30
nRCP.reps <- 30
nRCP.vector <- rep(c(nRCP.lower:nRCP.upper), nRCP.reps)


# run rcp models -----------------------------------------------------
for(i in 1:length(nRCP.vector)){
  
  nRCP <- nRCP.vector[i]
  run <- i
  curr.fit <- fit_rcp(complete_subset_list = complete_subset_list, nRCP = nRCP)
  fileName <- paste0("output/RegimixStats_nRCP", nRCP,"_run", run, ".RData")
  saveRDS(curr.fit, file = fileName)
  
}



