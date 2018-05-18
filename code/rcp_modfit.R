# RCP model fitting
#######


#---------------------------------------------------------#
# fit rcp models

make_rcp_dataobjs <- function(complete_subset_list){
  
  # load OTU table and covariates
  complete_subset_list <- readRDS(file = "derived_data/complete_subset_list.RData")
  covariates.species <- data.frame(complete_subset_list[["covariates"]], complete_subset_list[["otus"]])
  
  # define where abundance data starts in covariates.species, test with names(covariates.species)[n.abund]
  n.abund <- 17
  
  # choose species to model with
  # create a list of species with occurance > n to use for modelling
  species.n <- 20 # change this to desired minimum occurance
  species.count <- data.frame(count = sort(colSums(ifelse(covariates.species[,n.abund:ncol(covariates.species)]>0, 1, 0)), decreasing=T))
  species.count$species <- row.names(species.count)
  model.species.vector <- species.count$species[species.count$count>species.n]
  model.species.string <- paste0("cbind(", paste(model.species.vector, collapse=","),")")
  
  # choose which covariates to use
  model.covariates.vector <- names(covariates.species)[6:16] #CAREFUL!!
  model.covariates.string <- character(length(model.covariates.vector)*2)
  for (i in 1:length(model.covariates.vector)){
    pos2 <- i*2
    pos1 <- pos2-1
    model.covariates.string[pos1] <- paste0(model.covariates.vector[i],".1")
    model.covariates.string[pos2] <- paste0(model.covariates.vector[i],".2")
  }
  model.covariates.string <- paste0(model.covariates.string, collapse = "+")
  model.covariates.string
  
  #add size as a covariate
  model.covariates.vector <- c(names(covariates.species)[5], model.covariates.vector)
  model.covariates.string <- paste0(model.covariates.vector[1],"+", model.covariates.string)
  
  # calculate quadratic polynomial cols
  # you could do this via a lapply call for the poly() calcs, and use model.covariate.vector with some tricky renaming,
  # but i like to see the varaibles laid out for some reason (one check on the actual variables included is nice?).
  covar.data <- covariates.species[,model.covariates.vector]
  covar.data.forpredict <- covar.data
  covar.data <- data.frame(
    covar.data$size,
    poly(covar.data$density, 2),
    poly(covar.data$barkthick, 2),
    poly(covar.data$waterperc, 2),
    poly(covar.data$P, 2),
    poly(covar.data$K, 2),
    poly(covar.data$Ca, 2),
    poly(covar.data$Mn, 2),
    poly(covar.data$Fe, 2),
    poly(covar.data$Zn, 2),
    poly(covar.data$N, 2),
    poly(covar.data$C, 2))
  names(covar.data) <- c(
    "size",
    "density.1","density.2",
    "barkthick.1","barkthick.2",
    "waterperc.1","waterperc.2",
    "P.1","P.2",
    "K.1","K.2",
    "Ca.1","Ca.2",
    "Mn.1","Mn.2",
    "Fe.1","Fe.2",
    "Zn.1","Zn.2",
    "N.1","N.2",
    "C.1","C.2")
  ## convert categorical variables to factors if they exist
  covar.data$size <- factor(covar.data$size)
  
  # generate model data
  model.data <- data.frame(covariates.species[,model.species.vector], covar.data)
  
  # define model form
  RCP.form <- paste0(model.species.string,"~","1","+",model.covariates.string)
  # record site order
  site.names <- covariates.species$seq_sampName
  
  rcp.objs <- list(model.data = model.data, 
                   RCP.form = RCP.form, 
                   site.names = site.names,
                   model.covariates.vector = model.covariates.vector,
                   model.species.vector =  model.species.vector,
                   species.n = species.n,
                   covar.data.forpredict = covar.data.forpredict)
  return(rcp.objs)
  
}

fit_rcp <- function(complete_subset_list, nRCP){
  
  require(RCPmod)

  # load data
  rcp.objs <- make_rcp_dataobjs(complete_subset_list = complete_subset_list)
  
  # set control params
  my.cont <- list(maxit=2000, penalty=0.0001, penalty.tau=10, penalty.gamma=10)
  
  # run model
  tic <- proc.time()
  fit.regi <- regimix(form.RCP = rcp.objs$RCP.form, 
                      form.spp = NULL, 
                      data = rcp.objs$model.data, 
                      nRCP = nRCP, 
                      dist = "NegBin", 
                      control = my.cont, 
                      inits = "noPreClust", 
                      titbits = TRUE)
  toc <- proc.time()
  
  # write model fit stats
  modelStats <- list(sites = rcp.objs$site.names, 
                     covariates = rcp.objs$model.covariates.vector, 
                     species = rcp.objs$model.species.vector,
                     SppMin = rcp.objs$species.n, 
                     SppN = fit.regi$S, 
                     nRCP = fit.regi$nRCP, 
                     runtime = round((toc-tic)[3]/60),
                     AIC = fit.regi$AIC, 
                     BIC = fit.regi$BIC, 
                     postProbs = fit.regi$postProbs, 
                     logl = fit.regi$logl, 
                     coefs = fit.regi$coefs, 
                     penalties = unlist(my.cont), 
                     conv = fit.regi$conv)
  
  return(modelStats)
  
}

batch_fit_rcp <- function(complete_subset_list, nRCP.lower, nRCP.upper, nRCP.reps){
  
  #if you'd like to load data from a saved intermediate version of complete_subset_list instead...
  #complete_subset_list <- readRDS(file = "derived_data/complete_subset_list.Rdata")
  
  #specify a vector of nRCP values
  nRCP.vector <- rep(c(nRCP.lower:nRCP.upper), nRCP.reps)
  
  #for each run, fit the rcp model and save the output
  curr.fit.list <- list()
  for(i in 1:length(nRCP.vector)){
    
    curr.fit.list[[i]] <- fit_rcp(complete_subset_list = complete_subset_list, 
                                  nRCP = nRCP.vector[i])
    
  }
  names(curr.fit.list) <- paste('rcp', nRCP.vector, sep = "")
  
  #summarize results
  rcp_output_list <- curr.fit.list
  rcp_output_results <- data.frame(BIC = unlist(lapply(X = rcp_output_list, FUN = function(x){x$BIC})),
                                   AIC = unlist(lapply(X = rcp_output_list, FUN = function(x){x$AIC})),
                                   logLik = unlist(lapply(X = rcp_output_list, FUN = function(x){x$logl})),
                                   conv = unlist(lapply(X = rcp_output_list, FUN = function(x){x$conv})),
                                   minpp = unlist(lapply(X = rcp_output_list, FUN = function(x){min(colSums(x$postProbs))})),
                                   maxpp = unlist(lapply(X = rcp_output_list, FUN = function(x){max(colSums(x$postProbs))})),
                                   nRCP = unlist(lapply(X = rcp_output_list, FUN = function(x){x$nRCP})))
  
  #save df
  saveRDS(rcp_output_results, file = "derived_data/rcpFits/rcp_output_results_test.RData")
  
}


#---------------------------------------------------------#
# process cluster output

create_rcp_output_results <- function(){
  
  path <- "derived_data/rcpFits/fromCluster/"
  files <- list.files(path)
  
  rcpIntermed <- list()
  for(i in 1:length(files)){
    rcpIntermed[[i]] <- readRDS(file = paste0(path, files[i])) # each .RData file from the cluster is a run with fit.stats, collect these in list
  }
  
  #summarize results
  rcp_output_list <- rcpIntermed
  rcp_output_results <- data.frame(BIC = unlist(lapply(X = rcp_output_list, FUN = function(x){x$BIC})),
                                   AIC = unlist(lapply(X = rcp_output_list, FUN = function(x){x$AIC})),
                                   logLik = unlist(lapply(X = rcp_output_list, FUN = function(x){x$logl})),
                                   conv = unlist(lapply(X = rcp_output_list, FUN = function(x){x$conv})),
                                   minpp = unlist(lapply(X = rcp_output_list, FUN = function(x){min(colSums(x$postProbs))})),
                                   maxpp = unlist(lapply(X = rcp_output_list, FUN = function(x){max(colSums(x$postProbs))})),
                                   nRCP = unlist(lapply(X = rcp_output_list, FUN = function(x){x$nRCP})))
  
  #save df
  saveRDS(rcp_output_results, file = "derived_data/rcpFits/rcp_output_results.RData")
  
}

rcp_pick <- function(nRCP){
  
  #load data
  complete_subset_list <- readRDS(file = "derived_data/complete_subset_list.RData")
  rcp.objs <- make_rcp_dataobjs(complete_subset_list)
  
  #load fit params
  path <- "derived_data/rcpFits/fromCluster/"
  files <- list.files(path)
  pickedFiles <- files[grepl(paste0("nRCP",nRCP), files)] #there should be X reps with the same nRCP
  modelStats <- readRDS(file = paste0(path, pickedFiles[1])) #just grab the first one in the list
  
  # re-fit model w/o optimization
  modelStats$conv == 0
  if (modelStats$conv == 0) {
    my.cont <- list(penalty=0.0001, penalty.tau=10, penalty.gamma=10, optimise=FALSE) # no need to optimise, already have param values
    nRCP <- modelStats$nRCP
    params <- unlist(modelStats$coefs)
    model.data <- rcp.objs$model.data
    RCP.form <- rcp.objs$RCP.form
    fm_rcp <- regimix(form.RCP=RCP.form, form.spp=NULL, data=model.data, nRCP=nRCP,
                      dist="NegBin", control=my.cont, inits=params, titbits=TRUE)
    saveRDS(fm_rcp, file="derived_data/rcpFits/fm_rcp.RData")
    
  } else {
    message("model not converged, try again")
  }
  
}


#---------------------------------------------------------#
# load intermediates

load_rcp_output_results <- function(test){
  
  if(test == TRUE){
    rcp_output_results <- readRDS(file = "derived_data/rcpFits/rcp_output_results_test.RData")
  }else{
    rcp_output_results <- readRDS(file = "derived_data/rcpFits/rcp_output_results.RData")
  }
  
  return(rcp_output_results)
  
}

load_fm_rcp <- function(test){
  
  if(test == TRUE){
    fm_rcp <- readRDS(file = "derived_data/rcpFits/fm_rcp_test.RData")
  }else{
    fm_rcp <- readRDS(file = "derived_data/rcpFits/fm_rcp.RData")
  }
  
  return(fm_rcp)
}
