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
comm.otu.tmp <- load_matotu()
seqSamples <- load_seqSamples(mat.otu = comm.otu.tmp)
taxAndFunguild <- load_TaxAndFunguild(comm.otu.tmp = comm.otu.tmp)
comm.otu <- clean_comm(comm.otu.tmp = comm.otu.tmp, taxAndFunguild = taxAndFunguild)

# in-text results about the classification of OTUs...
# dim(taxAndFunguild)
# taxAndFunguild %>%
#   filter(genus != "unclassified") -> tmp
# num.genus.class <- dim(tmp)[1]
# tmp %>%
#   filter(Trophic.Mode != "unclassified") -> tmp2
# num.genus.class.guild<- dim(tmp2)[1]
# (num.genus.class.guild / num.genus.class) * 100

#saveRDS(taxAndFunguild, file = "derived_data/taxAndFunguild.RData") # derived_data/taxAndFunguild.RData
zanneTree <- load_zanne_tree()
traits.code <- mergeTraitData()  
complete_subset_list <- subset_and_make_list(seqSamples, traits.code, mat.otu = comm.otu)
save_complete_subset_list(complete_subset_list) # derived_data/complete_subset_list.RData
#dim(complete_subset_list$otus.trimmed)
#dim(complete_subset_list$otus)

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

# ** For testing **
# batch_fit_boral()

# ** For manuscript **
# commands associated with boral model fitting take > 24 hrs and were run separately with cluster computing
# script used for cluster computing: forCluster/clusterScript_woodEndophytes_boral.R
# output from cluster computing consists of one .RData file for each model fit
# to reproduce manuscript results...
# (1) download the zipped output files from this Google Drive link: 
# (2) place unzipped .RData files in this folder: derived_data/boralFits/fromCluster
# (3) extract summary data from .RData objects in the "fromCluster" folder to create smaller-sized .RData objects in the "boralFits" folder
#.... do so by uncommenting this section
#test.select <- FALSE
#extract_TraitLVs_allXs.Xcoefs.df(test.select) # makes Traits-and-LVs_allXs_allruns_Xcoefsreal.RData
#extract_TraitLVs_allXs.cor.df(test.select) # makes Traits-and-LVs_allXs_allruns_cordfreal.RData
#extract_TraitLVs_selectXs.cor.df(test.select) # makes Traits-and-LVs_selectXs_allruns_cordfreal.RData
#.... end of section

# you can also examine the mcmc chains, but do not upload this intermediate to GitHub because its too big (>2GB)
#### careful #### extract_TraitLVs_allXs.mcmc.obj(test.select) 
#TraitLVs_allXs.mcmc.obj <- readRDS("derived_data/boralFits/Traits-and-LVs_allXs_allruns_mcmcObjreal.RData")
#geweke.TraitsLVs_allX <- create_geweke_df(fit.list = TraitLVs_allXs.mcmc.obj, complete_subset_list, taxAndFunguild)
#summarizeGeweke_byParamType(geweke.df.ann = geweke.TraitsLVs_allX, modelLVonly = FALSE, allXs = TRUE) # output/boral_diagnostics/summarizeGeweke_byParamType_Traits-and-LVs.pdf = Fig S2

# (4) load minimal output from cluster computing...
TraitLVs_allXs.Xcoefs.df <- readRDS("derived_data/boralFits/Traits-and-LVs_allXs_allruns_Xcoefsreal.RData")
TraitLVs_allXs.cor.df <- readRDS("derived_data/boralFits/Traits-and-LVs_allXs_allruns_cordfreal.RData")
TraitLVs_selectXs.cor.df <- readRDS("derived_data/boralFits/Traits-and-LVs_selectXs_allruns_cordfreal.RData")

# # in-text results about number of OTU pairs used in analyses...
# fit.list <- TraitLVs_allXs.cor.df
# cor.list <- lapply(fit.list, function(x) x$cor.df)
# cor.df <- list_to_df(cor.list)
# cor.df.ann <- annotate_cor_withOTUInfo(cor.df = cor.df, taxAndFunguild = taxAndFunguild)
# cor.df.ann %>%
#   filter(source == "run10") -> cor.df.ann.curr
# cor.df.ann.curr %>%
#   mutate(otupair = paste0(otu1, otu2)) -> cor.df.ann.curr
# length(unique(cor.df.ann.curr$otupair)) == dim(cor.df.ann.curr)[1] # if TRUE, then each row is a unique OTU pair
# # total number of OTU pairs
# totalpairs <- dim(cor.df.ann.curr)[1]
# totalpairs
# cor.df.ann.curr %>%
#   group_by(pair_haveGenusIDs) %>%
#   summarize(n = length(otupair))
# # OTU pairs where both OTUs id'd to genus
# cor.df.ann.curr %>%
#   group_by(pair_samePhylum) %>%
#   summarize(n = length(otupair))
# # OTU pairs where both OTUs id'd to Asco or Basidio
# cor.df.ann.curr %>%
#   group_by(pair_sameTroph) %>%
#   summarize(n = length(otupair))
# # OTU pairs where both OTUs id'd to Sapro or Patho


# -------------------------------------------------------------------#
# boral_roleOf_traits_on_fungi.R : Summarize the direction and magnitude that wood traits explain OTU abundances
write_summary_Xcoefs(fit.list = TraitLVs_allXs.Xcoefs.df, allXs = TRUE) # output/boral_roleOfTraits/summary_Xcoefs.csv = Table 2
plot_summary_Xcoefs_byOTUId(fit.list = TraitLVs_allXs.Xcoefs.df, taxAndFunguild, allXs = TRUE) # output/boral_roleOfTraits/signifXcoefs_byOTUId_allX.pdf = Fig S6
# note for plot_summary_Xcoefs_byOTUId: OTU coef included if signif in 9 of 12 model runs or more

# -------------------------------------------------------------------#
# boral_cooccur.R : Summarize the direction and magnitude that OTU environmental and LV predictors covary
plot_cor_distributions(fit.list = TraitLVs_allXs.cor.df, taxAndFunguild, allXs = TRUE) # output/boral_cooccur/cor_distributions_allX.pdf = Fig 2
make_chordDiagrams(fit.list = TraitLVs_allXs.cor.df, taxAndFunguild, allXs = TRUE) # output/boral_cooccur/chordDiagrams.pdf = Fig 2
plot_cor_distributions(fit.list = TraitLVs_selectXs.cor.df, taxAndFunguild, allXs = FALSE) # output/boral_cooccur/cor_distributions_selectX.pdf = Fig S7
make_chordDiagrams_table(fit.list = TraitLVs_allXs.cor.df, taxAndFunguild, complete_subset_list, allXs = TRUE) # output/boral_cooccur/chordDiagrams.csv = Table S2
# note if the following functions are used for new data -- the y axes are fixed
plot_corFreq_phylo(fit.list = TraitLVs_allXs.cor.df, taxAndFunguild, allXs = TRUE) # output/boral_cooccur/corFreq_phylo.pdf = Fig 3
plot_corFreq_troph(fit.list = TraitLVs_allXs.cor.df, taxAndFunguild, allXs = TRUE) # output/boral_cooccur/corFreq_troph.pdf = Fig 3
#plot_cor_distributions_allruns(fit.list = TraitLVs_allXs.cor.df, taxAndFunguild, allXs = TRUE)

# in-text results about consistency of OTU pair correlations across model runs ...
# enviro.cor.btwRuns(fit.list = TraitLVs_allXs.cor.df, taxAndFunguild, allXs = TRUE) # %
# residual.cor.btwRuns(fit.list = TraitLVs_allXs.cor.df, taxAndFunguild, allXs = TRUE) # %

# in-text results that highlight notable shared enviro OTU correlations ...
enviro.out <- investigate_enviro_chordDiagram(fit.list = TraitLVs_allXs.cor.df, 
                                       fit.list2 = TraitLVs_allXs.Xcoefs.df, 
                                       taxAndFunguild, allXs = TRUE)
enviro.out$pair.rank
enviro.out$pair.info$ITSall_OTUa_5566ITSall_OTUd_35
taxAndFunguild %>%
  filter(OTUId == "ITSall_OTUd_35")

# in-text results that highlight notable residual OTU correlations ...
resid.out <- investigate_resid_chordDiagram(fit.list = TraitLVs_allXs.cor.df, 
                                            fit.list2 = TraitLVs_allXs.Xcoefs.df, 
                                            taxAndFunguild, allXs, 
                                            complete_subset_list, seqSamples, zanneTree)
resid.out$pair.rank
resid.out$pair.info[[8]]
taxAndFunguild %>%
  filter(OTUId == "ITSall_OTUb_7634")

taxAndFunguild %>%
  filter(grepl("trabicola", species))

# -------------------------------------------------------------------#
# rcp_modfit.R : Run mixtures-of-experts models to identify regions of common profiles (RCPs) among OTUs using wood trait measures as covariates

# ** For testing **
# (1) Run array of nRCPs and summarize stats and save as 'rcp_output_results'
# batch_fit_rcp(complete_subset_list = complete_subset_list, 
#               nRCP.lower = 2, nRCP.upper = 4, nRCP.reps = 1) # makes derived_data/rcpFits/rcp_output_results_test.RData
# (2) Re-run the 'best nRCP choice' and save the full details
# fm_rcp <- fit_rcp(complete_subset_list = complete_subset_list, 
#                   nRCP = 5)
# saveRDS(fm_rcp, file = "derived_data/rcpFits/fm_rcp_test.RData")

# ** For manuscript **
# commands associated with rcp model fitting take > 12 hrs and were run separately with cluster computing
# script used for cluster computing: forCluster/clusterScript_woodEndophytes_rcp.R
# to reproduce manuscript results...
# (1) download the zipped output files from this Google Drive link: 
# (2) place unzipped .RData files in this folder: derived_data/rcpFits/fromCluster
# (3) extract summary data from .RData objects in the "fromCluster" folder to create smaller-sized .RData objects in the "boralFits" folder
#.... do so by uncommenting this section
#create_rcp_output_results() # makes derived_data/rcpFits/rcp_output_results.RData
#.... end of section
# (4) save the .RData object that corresponds to the 'best nRCP choice'
#.... do so by uncommenting this section
#rcp_pick(nRCP = 5) # makes derived_data/rcpFits/fm_rcp.RData
#.... end of section

# (5) load minimal output from cluster computing...
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


