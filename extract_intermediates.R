#load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(GGally)
library(boral) #NOTE---- boral uses JAGS so you need to first install from http://sourceforge.net/projects/mcmc-jags/files/JAGS/

#load functions
source('code/boral_modfit.R')
source('code/boral_diagnostics.R')
source('code/load_data.R')

# extract minimal intermediates from boral fit objects and save as .RData files
test.select <- FALSE
#extract_TraitLVs_allXs.mcmc.obj(test.select) # don't upload this intermediate because its too big (>2GB)... see below
extract_TraitLVs_allXs.Xcoefs.df(test.select)
extract_TraitLVs_allXs.cor.df(test.select)
extract_TraitLVs_selectXs.cor.df(test.select)


# boral_diagnostics.R : Report fit diagnostics
TraitLVs_allXs.mcmc.obj <- readRDS("derived_data/boralFits/Traits-and-LVs_allXs_allruns_mcmcObjreal.RData")
complete_subset_list <- readRDS(file = "derived_data/complete_subset_list.RData")
taxAndFunguild <- readRDS(file = "derived_data/taxAndFunguild.RData")

# geweke convergence
geweke.TraitsLVs_allX <- create_geweke_df(fit.list = TraitLVs_allXs.mcmc.obj, complete_subset_list, taxAndFunguild)
summarizeGeweke_byParamType(geweke.df.ann = geweke.TraitsLVs_allX, modelLVonly = FALSE, allXs = TRUE) # output/boral_diagnostics/summarizeGeweke_byParamType_Traits-and-LVs.pdf = Fig S2

