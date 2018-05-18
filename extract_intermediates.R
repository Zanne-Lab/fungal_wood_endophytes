
source('code/boral_modfit.R')

# extract minimal intermediates from boral fit objects and save as .RData files
test.select <- FALSE
extract_TraitLVs_allXs.mcmc.obj(test.select)
extract_TraitLVs_allXs.Xcoefs.df(test.select)
extract_TraitLVs_allXs.cor.df(test.select)
extract_TraitLVs_selectXs.cor.df(test.select)
