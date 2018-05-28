#Linear version of remake.yml

#load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(GGally)
library(taxonlookup)
library(phytools)
library(diversitree)
library(vegan)
library(ggtree)
library(viridis)
library(mvabund)
library(boral) #NOTE---- boral uses JAGS so you need to first install from http://sourceforge.net/projects/mcmc-jags/files/JAGS/
library(circlize)
library(RCPmod)

#load functions
source("code/load_data.R")
source("code/plottingTheme.R")
source("code/summarize_OTUmat.R")
source("code/distances.R")
source("code/roleOf_traits_on_fungi.R")
source("code/boral_modfit.R")
source("code/boral_diagnostics.R")
source("code/boral_roleOf_traits_on_fungi.R")
source("code/boral_cooccur.R")
source("code/rcp_modfit.R")
source("code/rcp_diagnostics.R")
source("code/rcp_results.R")


# -------------------------------------------------------------------#
# load_data.R and plottingTheme.R : Load data and plotting theme
comm.otu <- load_matotu()
seqSamples <- load_seqSamples(mat.otu = comm.otu)
taxAndFunguild <- load_TaxAndFunguild(comm.otu)
#saveRDS(taxAndFunguild, file = "derived_data/taxAndFunguild.RData") # derived_data/taxAndFunguild.RData
zanneTree <- load_zanne_tree()
traits.code <- mergeTraitData()  
complete_subset_list <- subset_and_make_list(seqSamples, traits.code, mat.otu = comm.otu)
save_complete_subset_list(complete_subset_list) # derived_data/complete_subset_list.RData

# -------------------------------------------------------------------#
# summarize_OTUmat.R : Summarize the OTU matrix
calc_otu_summStats(comm.otu, taxAndFunguild) # output/otuSummary/otuSummary.csv
plot_sampleEffortCurves(comm.otu) # output/otuSummary/sampleEffortCurves.pdf = Fig S4
plot_richness(comm.otu, seqSamples, taxAndFunguild, zanneTree) # output/otuSummary/richness.pdf = Fig S1

# -------------------------------------------------------------------#
# distances.R : Compare wood phylogenetic distances, wood trait distances, and fungal community distances
pairDist <- calc_pairdist(zanneTree, comm.otu, traits.code, seqSamples)
plot_woodDistfungDist_bothContinuous(pairDist) # output/distances/fig_woodPhylo_woodTrait_fungDist.pdf = Fig 6
write_woodDistfungDist_summary(pairDist) # output/distances/summary_woodPhylo_woodTrait_fungDist.csv = Table S3

# -------------------------------------------------------------------#
# roleOf_traits_on_fungi.R : Find out which wood traits are most important in shaping the fungal community

## distance-based methods
plot_envfit_on_unconstrainedOrd(complete_subset_list, seqSamples) # output/roleOfTraits/unconstrained_envfit.pdf = Fig S5
summary_bestmodel_constrainedOrd(complete_subset_list) # output/roleOfTraits/constrained_bestmodelsummary.csv = Table S1

## mvabund methods
mod.m.list <- fit_mvabund(complete_subset_list)
write_mvabund_aicTable(mod.m.list) # output/roleOfTraits/mvabundAICs.csv = Table 1

# -------------------------------------------------------------------#
# boral_modfit.R : Run model-based constrained and unconstrained ordinations using the boral package

# input data = derived_data/complete_subset_list.RData

# ** For testing **
# batch_fit_boral() # only uncomment this if the input data changes
# output = derived_data/boralFits/test/*test.RData

# ** For manuscript **
# these processes take > 24 hrs and were run separately with cluster computing
# cluster scripts = ...
# output = derived_data/boralFits/fromCluster/*real.RData

# to load all output, unzip data above and uncomment the following functions
# test.select <- FALSE # change to false for testing
# LVonly <- load_boralIntermediates_lv(test = test.select)
# TraitsLVs_allXs <- load_boralIntermediates_lvenv(test = test.select, allXs = TRUE)
# TraitsLVs_selectXs <- load_boralIntermediates_lvenv(test = test.select, allXs = FALSE)

# to load minimal output...
TraitLVs_allXs.Xcoefs.df <- readRDS("derived_data/boralFits/Traits-and-LVs_allXs_allruns_Xcoefsreal.RData")
TraitLVs_allXs.cor.df <- readRDS("derived_data/boralFits/Traits-and-LVs_allXs_allruns_cordfreal.RData")
TraitLVs_selectXs.cor.df <- readRDS("derived_data/boralFits/Traits-and-LVs_selectXs_allruns_cordfreal.RData")

# -------------------------------------------------------------------#
# boral_roleOf_traits_on_fungi.R : Summarize the direction and magnitude that wood traits explain OTU abundances
write_summary_Xcoefs(fit.list = TraitLVs_allXs.Xcoefs.df, allXs = TRUE) # output/boral_roleOfTraits/summary_Xcoefs.csv = Table 2
plot_summary_Xcoefs_byOTUId(fit.list = TraitLVs_allXs.Xcoefs.df, taxAndFunguild, allXs = TRUE) # output/boral_roleOfTraits/signifXcoefs_byOTUId_allX.pdf = Fig S6

# -------------------------------------------------------------------#
# boral_cooccur.R : Summarize the direction and magnitude that OTU environmental and LV predictors covary
plot_cor_distributions(fit.list = TraitLVs_allXs.cor.df, taxAndFunguild, allXs = TRUE) # output/boral_cooccur/cor_distributions_allX.pdf = Fig 2
make_chordDiagrams(fit.list = TraitLVs_allXs.cor.df, taxAndFunguild, allXs = TRUE) # output/boral_cooccur/chordDiagrams.pdf = Fig 2
plot_cor_distributions(fit.list = TraitLVs_selectXs.cor.df, taxAndFunguild, allXs = FALSE) # output/boral_cooccur/cor_distributions_selectX.pdf = Fig S7
make_chordDiagrams_table(fit.list = TraitLVs_allXs.cor.df, taxAndFunguild, complete_subset_list, allXs = TRUE) # output/boral_cooccur/chordDiagrams.csv = Table S2
plot_corFreq_phylo(fit.list = TraitLVs_allXs.cor.df, taxAndFunguild, allXs = TRUE) # output/boral_cooccur/corFreq_phylo.pdf = Fig 3
plot_corFreq_troph(fit.list = TraitLVs_allXs.cor.df, taxAndFunguild, allXs = TRUE) # output/boral_cooccur/corFreq_troph.pdf = Fig 3



# -------------------------------------------------------------------#
# rcp_modfit.R : Run mixtures-of-experts models to identify regions of common profiles (RCPs) among OTUs using wood trait measures as covariates

# ** For testing **
# 1. Run array of nRCPs and summarize stats and save as 'rcp_output_results'
# batch_fit_rcp(complete_subset_list = complete_subset_list, 
#               nRCP.lower = 2, nRCP.upper = 4, nRCP.reps = 1)
# output = derived_data/rcpFits/rcp_output_results_test.RData
# 2. Re-run the 'best nRCP choice' and save the full details
# fm_rcp <- fit_rcp(complete_subset_list = complete_subset_list, 
#                   nRCP = 5)
# saveRDS(fm_rcp, file = "derived_data/rcpFits/fm_rcp_test.RData")

# ** For manuscript **
# 1. Run array of nRCPs; these processes take a while and were run separately with cluster computing
# cluster scripts = ... clusterScript_rcp.R
# output = derived_data/rcpFits/fromCluster....RData
# 2. Summarize stats from all .RData objects and save 'rcp_output_results'
# create_rcp_output_results()
# output = derived_data/rcpFits/fromCluster/rcp_output_results.RData

# 3. Save the .RData object that corresponds to the 'best nRCP choice'
# rcp_pick(nRCP = 5) 
# output = derived_data/rcpFits/fromCluster/fm_rcp.RData

# Load the intermediate data into this workflow
rcp_output_results <- load_rcp_output_results(test = FALSE)
fm_rcp <- load_fm_rcp(test = FALSE)


# -------------------------------------------------------------------#
# rcp_diagnostics.R : Report fit diagnostics
plot_numOfgroups_diagnostics(rcp_output_results = rcp_output_results) # output/rcp_diagnostics/numOfgroups_diagnostics.pdf = Fig S3

# -------------------------------------------------------------------#
# rcp_results.R 
plot_confus(fm_rcp = fm_rcp, complete_subset_list = complete_subset_list, seqSamples = seqSamples, zanneTree = zanneTree) # output/rcp_results/species_confmat.pdf = Fig 5
plot_studyDesign(seqSamples = seqSamples, zanneTree = zanneTree) # output/rcp_results/species_studydesign.pdf = Fig 5
plot_partialEffects_select4(complete_subset_list = complete_subset_list, fm_rcp = fm_rcp) # output/rcp_results/partialEffects_select4.pdf = Fig 4
plot_partialEffect_theRest1(complete_subset_list = complete_subset_list, fm_rcp = fm_rcp) # output/rcp_results/partialEffects_all.pdf = Fig S9
plot_partialEffect_theRest2(complete_subset_list = complete_subset_list, fm_rcp = fm_rcp) # output/rcp_results/partialEffects_all.pdf = Fig S9
plot_mostabundTaxa_perRCP(fm_rcp, taxAndFunguild) # output/rcp_results/expectedAbund_perRCP.pdf = Fig S8


