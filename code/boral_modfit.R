# boral - model fitting
#NOTE---- boral uses JAGS so you need to first install JAGS-3.0.0.exe or higher from http://sourceforge.net/projects/mcmc-jags/files/JAGS/
#######

# mean model:
# for sample i and OTU j ...
# g(μ_{ij}) = alpha_i + beta_{0j} + \bm{x}^T_i\bm{beta}_j + \bm{z}^T_i\bm{theta}_j
# where g() is a log-link

# model components:
# alpha_i = random site effects = row.coefs.ID1 (samples)
# beta_{0j} = OTU-specific intercepts = lv.coefs.1
# \bm{x}^T_i = X covariate data
# \bm{beta}_j = OTU-specific X coefficients = X.coefs
# \bm{z}^T_i = Latent variables = lvs.1 and lvs.2
# \bm{theta}_j = OTU-specific LV coefficients = lv.coefs.2 and lv.coefs.3

# multiplicative random effect:
# u[i,j] ~ dgamma(1/lv.coefs[j,num.lv+2], 1/lv.coefs[j,num.lv+2])
# lv.coefs[j,num.lv+2] = Dispersion = lv.coefs.4

#---------------------------------------------------------#
# set up

make_boral_dataobjs_allX <- function(complete_subset_list){
  
  #identify the OTU matrix and covariates
  mat.otu<-complete_subset_list[["otus.trimmed"]]
  xVars<-complete_subset_list[["covariates"]]
  
  #identify the environmental covariates
  xVars %>%
    select(size, density, barkthick, waterperc,
           P, K, Ca, Mn, Fe, Zn, N, C) -> woodTraits
  
  #mean center continuous variables
  select <- colnames(woodTraits) != "size"
  woodTraits.contin <- scale(woodTraits[,select], scale = FALSE)
  
  #recode any categorical predictors as a dummy variable
  dummy <- data.frame(model.matrix( ~ size - 1, data = woodTraits))
  woodTraits %>%
    mutate(sizesmall = dummy$sizesmall) %>%
    select(-size) -> woodTraits
  
  #put the categorical and continous variables together
  woodTraits.mc <- data.frame(sizesmall = woodTraits$sizesmall, woodTraits.contin)
  
  #no groups of row effects
  row.ids <- NA
  
  boral_dataobjs <- list(Y = mat.otu, X = woodTraits.mc, rowEffs = row.ids)
  return(boral_dataobjs)
}

make_boral_dataobjs_selectX <- function(complete_subset_list){
  
  #identify the OTU matrix and covariates
  mat.otu<-complete_subset_list[["otus.trimmed"]]
  xVars<-complete_subset_list[["covariates"]]
  
  #identify the wood trait covariates
  #decided not to include density; it is highly correlated with waterperc
  xVars %>%
    select(size, C, waterperc, barkthick) -> woodTraits
  
  #mean center continuous variables
  select <- colnames(woodTraits) != "size"
  woodTraits.contin <- scale(woodTraits[,select], scale = FALSE)
  
  #recode any categorical predictors as a dummy variable
  dummy <- data.frame(model.matrix( ~ size - 1, data = woodTraits))
  woodTraits %>%
    mutate(sizesmall = dummy$sizesmall) %>%
    select(-size) -> woodTraits
  
  #put the categorical and continous variables together
  woodTraits.mc <- data.frame(sizesmall = woodTraits$sizesmall, woodTraits.contin)
  
  #no groups of row effects
  row.ids <- NA
  
  boral_dataobjs <- list(Y = mat.otu, X = woodTraits.mc, rowEffs = row.ids)
  return(boral_dataobjs)
}

set_mcmc_controls <- function(seed, runType){
  
  if(runType == "test"){
    n.burnin = 50
    n.iteration = 100
    n.thin = 1
  }
  
  if(runType == "real"){
    n.burnin = 100000 
    n.iteration = 200000
    n.thin = 10
  }
  
  mcmc.controls <- list(n.burnin = n.burnin,
       n.iteration = n.iteration,
       n.thin = n.thin,
       seed = seed)
  
  return(mcmc.controls)
}

prior_controls <- function(){
  
  prior.controls = list(type = c("normal","normal","normal","uniform"),
                        hypparams = c(100, #OTU-specific intercepts and row effects
                                      20,  #lv coefs
                                      100, #X coefs
                                      20)) #dispersion
}


#---------------------------------------------------------#
# extract info from the boral object

extract_lvsplotCoords <- function(fit, Y){
  
  lvsplot.vals <- lvsplot(fit, biplot = TRUE, est = "median", return.vals = TRUE)
  scaled.lvs <- data.frame(seq_sampName = row.names(Y),
                           LV1 = lvsplot.vals$scaled.lvs[,1], 
                           LV2 = lvsplot.vals$scaled.lvs[,2]) # LV coordinates for each sample
  scaled.lv.coefs <- data.frame(OTUId = colnames(Y),
                                LV1.coef = lvsplot.vals$scaled.lv.coefs[,1], 
                                LV2.coef = lvsplot.vals$scaled.lv.coefs[,2]) # LV coefs for each OTU, e.g. to make a biplot
  
  lvsplot.data <- list(scaled.lvs = scaled.lvs, scaled.lv.coefs = scaled.lv.coefs)
  
  return(lvsplot.data)
  
}

extract_residVfitted <- function(fit){
  
  # residuals
  ds.residuals <- data.frame(ds.residuals(fit, est = "median")$residuals)
  ds.residuals$seq_sampName <- row.names(ds.residuals)
  ds.residuals %>%
    gather(key = "OTUId", value = "residuals", -seq_sampName) -> residuals.df
  
  # fitted values
  fitted.vals <- data.frame(fitted(fit, est = "median")$out)
  fitted.vals$seq_sampName <- row.names(fitted.vals)
  fitted.vals %>%
    gather(key = "OTUId", value = "fitted.vals", -seq_sampName) -> fitted.vals.df
  
  # join them
  residuals.df %>%
    left_join(fitted.vals.df) -> resid.data
  
  return(resid.data)
  
}

extract_Xcoefs<-function(fit){
  
  #make dfs
  OTUId <- row.names(fit$X.coefs.median)
  X.df <- data.frame(OTUId = OTUId, fit$X.coefs.median)
  upper.df <- data.frame(OTUId = OTUId, fit$hpdintervals$X.coefs[,,"upper"])
  lower.df <- data.frame(OTUId = OTUId, fit$hpdintervals$X.coefs[,,"lower"])
  
  #make long
  X.df %>%
    gather(key = "term", value = "X", -OTUId) -> X.df.l
  upper.df %>%
    gather(key = "term", value = "X.upper", -OTUId) -> upper.df.l
  lower.df %>%
    gather(key = "term", value = "X.lower", -OTUId) -> lower.df.l
  
  #join dfs
  X.df.l %>%
    left_join(upper.df.l) %>%
    left_join(lower.df.l) -> Xcoefs.df
  
  return(Xcoefs.df)
  
}

extract_cors<-function(corobj, corType){
  
  # Full correlation matrix
  cor <- corobj$cor
  cor.df <- extract_uniquePairDists(cor) #make into a long df
  colnames(cor.df) <- c("otu1","otu2", "cor")
  
  # A correlation matrix containing only the “significant" correlations 
  # whose 95% highest posterior interval does not contain zero. 
  # All non-significant correlations set to zero.
  sig.cor <- corobj$sig.cor
  sig.cor.df <- extract_uniquePairDists(sig.cor) #make into a long df
  colnames(sig.cor.df) <- c("otu1","otu2", "sig.cor")
  
  # join dfs and id 'signif' cells
  cor.df %>%
    left_join(sig.cor.df) %>%
    mutate(signif = ifelse(sig.cor == 0,
                           "no", "yes")) %>%
    select(otu1, otu2, cor, signif) -> df
  
  # replace signif with more specific levels
  df[df$signif == "yes" & df$cor > 0, "signif"] <- "positive"
  df[df$signif == "yes" & df$cor < 0, "signif"] <- "negative"
  
  #update colnames to reflect corType
  colnames(df)[3:4] <- paste(corType, colnames(df)[3:4], sep="_")
  
  return(df)
  
}


#---------------------------------------------------------#
# fit functions

fit_boral_lvOnly <- function(complete_subset_list, seed, runType){
  
  require(boral)
  
  mcmc_controls <- set_mcmc_controls(seed = seed, runType = runType)
  boral_dataobjs <- make_boral_dataobjs_allX(complete_subset_list)  # uses "_allX"
  Y <- boral_dataobjs[["Y"]]
  prior.controls <- prior_controls()
  
  fit <- boral(y = Y,
               family = "negative.binomial",
               num.lv = 2,
               row.eff = "random",
               calc.ics = TRUE, # for model selection
               save.model = TRUE, # to access mcmc samples
               mcmc.control = mcmc_controls,
               prior.control = prior.controls)
  
  #extract and save objects from fit ...
  
  #extract the residual covariance matrix and just keep the trace
  rescor <- get.residual.cor(fit)
  
  #mcmc object
  mcmc.obj <- as.mcmc.list(fit$jags.model$BUGSoutput)
  
  #lvsplot coordinates
  lvsplot.data <- extract_lvsplotCoords(fit, Y)
  
  #residual vs fitted data
  resid.data <- extract_residVfitted(fit)
  #plot(fit)
  
  #row coef estimates
  rowcoef.mean<-list()
  for(i in 1:length(names(fit$row.coefs))){
    tmp <- fit$row.coefs[[i]]
    rowcoef.mean[[i]] <- tmp$mean
  }
  names(rowcoef.mean) <- names(fit$row.coefs)
  
  #put it all together
  fit.stats <- list(rescor.trace = rescor$trace,
                    mcmc.obj = mcmc.obj,
                    lvsplot.data = lvsplot.data, 
                    resid.data = resid.data,
                    lv.mean = fit$lv.mean,
                    lv.coefs.mean = fit$lv.coefs.mean,
                    rowcoef.mean = rowcoef.mean,
                    ics = fit$ics)
  
  return(fit.stats)
  
}

fit_boral_allX <- function(complete_subset_list, seed, runType){
  
  require(boral)
  
  mcmc_controls <- set_mcmc_controls(seed, runType = runType)
  boral_dataobjs <- make_boral_dataobjs_allX(complete_subset_list)
  Y <- boral_dataobjs[["Y"]]
  X <- boral_dataobjs[["X"]]
  prior.controls <- prior_controls()
  
  fit <- boral(y = Y,
               X = X,
               family = "negative.binomial",
               num.lv = 2,
               row.eff = "random",
               calc.ics = TRUE, # for model selection
               save.model = TRUE, # to access mcmc samples
               mcmc.control = mcmc_controls,
               prior.control = prior.controls)
  
  #extract and save objects from fit ...
  
  #mcmc object
  mcmc.obj <- as.mcmc.list(fit$jags.model$BUGSoutput)
  
  #lvsplot coordinates
  lvsplot.data <- extract_lvsplotCoords(fit, Y)
  
  #residual vs fitted data
  resid.data <- extract_residVfitted(fit)
  
  #row coef estimates
  rowcoef.mean <- list()
  for(i in 1:length(names(fit$row.coefs))){
    tmp <- fit$row.coefs[[i]]
    rowcoef.mean[[i]] <- tmp$mean
  }
  names(rowcoef.mean) <- names(fit$row.coefs)
  
  # .... data unique to modelID == "Traits-and-LVs"
  
  # X coefs
  Xcoefs.df <- extract_Xcoefs(fit)
  # determine if 0 is within the highest posterior density interval aka credible interval
  Xcoefs.df %>%
    mutate(signif = ifelse(X.upper < 0 & X.lower < 0 | X.upper > 0 & X.lower > 0,
                           "yes", "no")) -> Xcoefs.df
  
  # shared environment and residual correlations between OTUs
  enviro.cor <- get.enviro.cor(fit, prob = 0.95)
  enviro.cor.df <- extract_cors(corobj = enviro.cor, corType = "enviro")
  residual.cor <- get.residual.cor(fit, prob = 0.95)
  residual.cor.df <- extract_cors(corobj = residual.cor, corType = "residual")
  enviro.cor.df %>%
    left_join(residual.cor.df) -> cor.df
  
  #put it all together
  fit.stats <- list(rescor.trace = residual.cor$trace,
                    mcmc.obj = mcmc.obj,
                    lvsplot.data = lvsplot.data, 
                    resid.data = resid.data,
                    Xcoefs.df = Xcoefs.df,
                    cor.df = cor.df,
                    lv.mean = fit$lv.mean,
                    lv.coefs.mean = fit$lv.coefs.mean,
                    rowcoef.mean = rowcoef.mean,
                    ics = fit$ics)
  
  return(fit.stats)
  
}

fit_boral_selectX <- function(complete_subset_list, seed, runType){
  
  require(boral)
  
  mcmc_controls <- set_mcmc_controls(seed, runType = runType)
  boral_dataobjs <- make_boral_dataobjs_selectX(complete_subset_list)
  Y <- boral_dataobjs[["Y"]]
  X <- boral_dataobjs[["X"]]
  prior.controls <- prior_controls()
  
  fit <- boral(y = Y,
               X = X,
               family = "negative.binomial",
               num.lv = 2,
               row.eff = "random",
               calc.ics = TRUE, # for model selection
               save.model = TRUE, # to access mcmc samples
               mcmc.control = mcmc_controls,
               prior.control = prior.controls)
  
  #extract and save objects from fit ...
  
  #mcmc object
  mcmc.obj <- as.mcmc.list(fit$jags.model$BUGSoutput)
  
  #lvsplot coordinates
  lvsplot.data <- extract_lvsplotCoords(fit, Y)
  
  #residual vs fitted data
  resid.data <- extract_residVfitted(fit)
  
  #row coef estimates
  rowcoef.mean <- list()
  for(i in 1:length(names(fit$row.coefs))){
    tmp <- fit$row.coefs[[i]]
    rowcoef.mean[[i]] <- tmp$mean
  }
  names(rowcoef.mean) <- names(fit$row.coefs)
  
  # .... data unique to modelID == "Traits-and-LVs"
  
  # X coefs
  Xcoefs.df <- extract_Xcoefs(fit)
  # determine if 0 is within the highest posterior density interval aka credible interval
  Xcoefs.df %>%
    mutate(signif = ifelse(X.upper < 0 & X.lower < 0 | X.upper > 0 & X.lower > 0,
                           "yes", "no")) -> Xcoefs.df
  
  # shared environment and residual correlations between OTUs
  enviro.cor <- get.enviro.cor(fit, prob = 0.95)
  enviro.cor.df <- extract_cors(corobj = enviro.cor, corType = "enviro")
  residual.cor <- get.residual.cor(fit, prob = 0.95)
  residual.cor.df <- extract_cors(corobj = residual.cor, corType = "residual")
  enviro.cor.df %>%
    left_join(residual.cor.df) -> cor.df
  
  #put it all together
  fit.stats <- list(rescor.trace = residual.cor$trace,
                    mcmc.obj = mcmc.obj,
                    lvsplot.data = lvsplot.data, 
                    resid.data = resid.data,
                    Xcoefs.df = Xcoefs.df,
                    cor.df = cor.df,
                    lv.mean = fit$lv.mean,
                    lv.coefs.mean = fit$lv.coefs.mean,
                    rowcoef.mean = rowcoef.mean,
                    ics = fit$ics)
  
  return(fit.stats)
  
}

batch_fit_boral <- function(){
  
  # load data -----------------------------------------------------
  complete_subset_list <- readRDS(file = "derived_data/complete_subset_list.RData")
  runType <- "test"
  
  # LV-only model -----------------------------------------------------
  NUMRUNS <- 3
  SEED <- .Random.seed[1:NUMRUNS] #set seed for different runs
  for(i in 1:NUMRUNS){
     curr.fit.stats <- fit_boral_lvOnly(complete_subset_list, seed = SEED[i], runType = runType)
     curr.fileName <- paste0("derived_data/boralFits/test/LV-only_run", i, "_runType", runType, ".RData")
     saveRDS(curr.fit.stats, file = curr.fileName)
  }
  
  # Traits-and-LVs model (all X vars) -----------------------------------------------------
  NUMRUNS <- 3
  SEED <- .Random.seed[1:NUMRUNS] #set seed for different runs
  for(i in 1:NUMRUNS){
    curr.fit.stats <- fit_boral_allX(complete_subset_list, seed = SEED[i], runType = runType)
    curr.fileName <- paste0("derived_data/boralFits/test/Traits-and-LVs_allX_run", i, "_runType", runType, ".RData")
    saveRDS(curr.fit.stats, file = curr.fileName)
  }
  
  # Traits-and-LVs model (select X vars) -----------------------------------------------------
  NUMRUNS <- 3
  SEED <- .Random.seed[1:NUMRUNS] #set seed for different runs
  for(i in 1:NUMRUNS){
    curr.fit.stats <- fit_boral_selectX(complete_subset_list, seed = SEED[i], runType = runType)
    curr.fileName <- paste0("derived_data/boralFits/test/Traits-and-LVs_selectX_run", i, "_runType", runType, ".RData")
    saveRDS(curr.fit.stats, file = curr.fileName)
  }
  
}


#---------------------------------------------------------#
# load boral fit.list objects

read.in.boralobjs <- function(path, file.names){
  
  fit.list <- list()
  for(i in 1:length(file.names)){
    fit.list[[i]] <- readRDS(file = paste0(path, file.names[i])) 
  }
  names(fit.list) <- paste0("run", 1:length(fit.list))
  
  return(fit.list)
  
}

load_boralIntermediates_lv <- function(test){
  
  #identify the folder where the files are located
  if(test == TRUE){
    path <- "derived_data/boralFits/test/"
  }else{
    path <- "derived_data/boralFits/fromCluster/"
  }
  
  #identify the files
  files <- list.files(path)
  file_run <- files[grepl("LV-only", files)]
  
  #read in the files and save them in a named list
  boralIntermed <- read.in.boralobjs(path = path, file.names = file_run)
  
  return(boralIntermed)
  
}

load_boralIntermediates_lvenv <- function(test, allXs){
  
  #identify the folder where the files are located
  if(test == TRUE){
    path <- "derived_data/boralFits/test/"
  }else{
    path <- "derived_data/boralFits/fromCluster/"
  }
  
  #identify the files
  files <- list.files(path)
  if(allXs == TRUE){
    file_run <- files[grepl("Traits-and-LVs_allX", files)]
    }else{
    file_run <- files[grepl("Traits-and-LVs_selectX", files)]
  }
  
  #read in the files and save them in a named list
  boralIntermed <- read.in.boralobjs(path = path, file.names = file_run)
  
  return(boralIntermed)
  
}


#---------------------------------------------------------#
# extract, save, load minimal intermediates from boral fit.list objects

extract_TraitLVs_allXs.mcmc.obj <- function(test.select){
  
  TraitsLVs_allXs <- load_boralIntermediates_lvenv(test = test.select, allXs = TRUE)
  TraitLVs_allXs.mcmc.obj <- lapply(TraitsLVs_allXs, function(x){
    list(mcmc.obj = x$mcmc.obj)
  })
  
  if(test.select == TRUE){
    runType <- "test"
  }else{
    runType <- "real"
  }
  
  curr.fileName <- paste0("derived_data/boralFits/Traits-and-LVs_allXs_allruns_mcmcObj", runType, ".RData")
  saveRDS(TraitLVs_allXs.mcmc.obj, file = curr.fileName)
}

extract_TraitLVs_allXs.Xcoefs.df <- function(test.select){
  
  TraitsLVs_allXs <- load_boralIntermediates_lvenv(test = test.select, allXs = TRUE)
  TraitLVs_allXs.Xcoefs.df <- lapply(TraitsLVs_allXs, function(x){
    list(Xcoefs.df = x$Xcoefs.df)
  })
  
  if(test.select == TRUE){
    runType <- "test"
  }else{
    runType <- "real"
  }
  
  curr.fileName <- paste0("derived_data/boralFits/Traits-and-LVs_allXs_allruns_Xcoefs", runType, ".RData")
  saveRDS(TraitLVs_allXs.Xcoefs.df, file = curr.fileName)

}

extract_TraitLVs_allXs.cor.df <- function(test.select){
  
  TraitsLVs_allXs <- load_boralIntermediates_lvenv(test = test.select, allXs = TRUE)
  TraitLVs_allXs.cor.df <- lapply(TraitsLVs_allXs, function(x){
    list(cor.df = x$cor.df)
  })
  
  if(test.select == TRUE){
    runType <- "test"
  }else{
    runType <- "real"
  }
  
  curr.fileName <- paste0("derived_data/boralFits/Traits-and-LVs_allXs_allruns_cordf", runType, ".RData")
  saveRDS(TraitLVs_allXs.cor.df, file = curr.fileName)
}

extract_TraitLVs_selectXs.cor.df <- function(test.select){
  
  TraitsLVs_selectXs <- load_boralIntermediates_lvenv(test = test.select, allXs = FALSE)
  TraitLVs_selectXs.cor.df <- lapply(TraitsLVs_selectXs, function(x){
    list(cor.df = x$cor.df)
  })
  
  if(test.select == TRUE){
    runType <- "test"
  }else{
    runType <- "real"
  }
  
  curr.fileName <- paste0("derived_data/boralFits/Traits-and-LVs_selectXs_allruns_cordf", runType, ".RData")
  saveRDS(TraitLVs_selectXs.cor.df, file = curr.fileName)
  
}




