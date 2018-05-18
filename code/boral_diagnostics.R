# boral - model fitting diagnostics
#######

#---------------------------------------------------------#
# summarize parameter convergence

create_geweke_df <- function(fit.list, complete_subset_list, taxAndFunguild){
  
  #Geweke (1992) proposed a convergence diagnostic for Markov chains based on a test for equality of the means of the first and last part of a Markov chain (by default the first 10% and the last 50%). 
  #If the samples are drawn from the stationary distribution of the chain, the two means are equal and Geweke's statistic has an asymptotically standard normal distribution
  #if the p-value is less than 0.05, there is evidence the MCMC chain has not converged
  
  # select the mcmc obj and calculate geweke
  geweke.list <- lapply(fit.list, function(x){
    
    mcmc.obj <- x$mcmc.obj
    geweke <- geweke.diag(mcmc.obj, frac1=0.1, frac2=0.5)
    geweke.z <- geweke[[1]]$z
    geweke.df <- data.frame(param = names(geweke.z), z = geweke.z)
    
    return(geweke.df)
  })
  
  # reformat, convert z to pval, do holm pval adjustment
  geweke.df <- list_to_df(geweke.list)
  geweke.df %>%
    mutate(pval = 2*pnorm(abs(z), lower.tail = FALSE)) %>%
    mutate(adj.pval = p.adjust(pval, method = "holm")) %>%
    mutate(converg = ifelse(adj.pval > 0.05, "yes", "no")) -> geweke.df
  
  # annotate params
  geweke.df %>%
    separate(param, into = c("paramType","string"), sep = "\\[", fill = "right", remove = F) %>%
    separate(string, into = c("paramNum1","string2"), sep = "\\,", fill = "right") %>%
    separate(string2, into = c("paramNum2", "drop"), sep = "\\]", fill = "right") %>%
    select(source, param, paramType, paramNum1, paramNum2, adj.pval, converg) -> geweke.df.ann
  
  #fix paramNum1 for rowIDs
  selection <- geweke.df.ann$paramType %in% c("row.coefs.ID1", "row.coefs.ID2")
  vec <-geweke.df.ann[selection,"paramNum1"]
  geweke.df.ann[selection, "paramNum1"] <- unlist(strsplit(vec, "\\]"))
  
  # param structure
  geweke.df.ann %>%
    group_by(paramType) %>%
    summarize(length(unique(paramNum1)),
              length(unique(paramNum2)))
  
  # add paramTypeNum category
  geweke.df.ann %>%
    mutate(paramTypeNum = paste0(paramType, paramNum2)) -> geweke.df.ann
  
  # annotate OTUs for lv.coefs and X.coefs
  Y <- complete_subset_list[['otus.trimmed']]
  OTUindx <- data.frame(OTUId = colnames(Y), 
                        OTUnum = as.character(1:length(colnames(Y))),
                        numSamps_present = apply(Y > 0, 2, sum))
  taxa.indx <- taxAndFunguild[,c("OTUId","OTUId_ann")]
  OTUindx %>%
    left_join(taxa.indx) -> OTUindx
  geweke.df.ann %>%
    filter(paramType %in% c("lv.coefs", "X.coefs")) %>%
    mutate(OTUnum = paramNum1) %>%
    left_join(OTUindx) -> geweke.lvcoefs.X
  
  # annotate sample ids for lvs and row.coefs.ID1
  covariates <- complete_subset_list[["covariates"]]
  sampindx <- data.frame(seq_sampName = covariates$seq_sampName, 
                        sampNum = as.character(1:dim(covariates)[1]))
  geweke.df.ann %>%
    filter(paramType %in% c("lvs","row.coefs.ID1")) %>%
    mutate(sampNum = paramNum1) %>%
    left_join(sampindx) -> geweke.lvs.rowID1
  
  # join dfs
  geweke.df.ann %>%
    filter(!paramType %in% c("lv.coefs","lvs","row.coefs.ID1", "X.coefs")) %>%
    full_join(geweke.lvcoefs.X) %>%
    full_join(geweke.lvs.rowID1) -> geweke.df.ann
  
  # remove the deviance
  geweke.df.ann %>%
    filter(paramType != "deviance") -> geweke.df.ann
  
  return(geweke.df.ann)
}

summarizeGeweke_byParamType <- function(geweke.df.ann, modelLVonly, allXs){
  
  if(modelLVonly == TRUE){
    modelID <- "LV-only"
  }else{
    modelID <- "Traits-and-LVs"
  }
  
  if(allXs == TRUE){
    modelID2 <- "allX"
    xcoef.order <- c("sizesmall","density","barkthick","waterperc","P","K","Ca","Mn","Fe","Zn","N","C")
  }else{
    modelID2 <- "selectX"
    xcoef.order <- c("sizesmall","C","waterperc","barkthick")
  }
  
  # remove lv.coefs[1,3]
  selection <- geweke.df.ann$paramType == "lv.coefs" & geweke.df.ann$paramNum1 == "1" & geweke.df.ann$paramNum2 == "3"
  geweke.df.ann <- geweke.df.ann[!selection,]

  # summarize by parameter type and model run
  geweke.df.ann %>%
    group_by(paramTypeNum, source) %>%
    summarize(total = length(converg),
      numConverged = sum(converg == "yes"),
      percConverged = (numConverged/total)*100) -> summ
  
  # average percent converged aggregated by run
  n.runs <- length(unique(summ$source))
  summ %>%
    group_by(paramTypeNum) %>%
    summarize(mean = mean(percConverged),
              sd = sd(percConverged),
              se = sd / sqrt(n.runs)) -> summ2
  
  # add pretty names
  summ2 %>%
    filter(grepl("X.coefs", paramTypeNum)) %>%
    separate(paramTypeNum, into = c("drop","num"), sep = "X.coefs") %>%
    mutate(num = as.numeric(num)) %>%
    select(-drop) %>%
    arrange(num) -> summ2.xcoefs
  summ2.xcoefs <- cbind(summ2.xcoefs, xcoef.order)
  summ2.xcoefs %>%
    mutate(paramType = "Xcoef") %>%
    mutate(paramTypeNum = paste0(paramType, "_", xcoef.order)) %>%
    select(paramTypeNum, mean, se) -> summ2.xcoefs
  summ2 %>%
    filter(!grepl("X.coefs", paramTypeNum)) %>%
    mutate(paramTypeNum = ifelse(paramTypeNum == "row.coefs.ID1NA", "row.coefs",paramTypeNum)) %>%
    mutate(paramTypeNum = ifelse(paramTypeNum == "row.sigma.ID1NA", "row.sigma",paramTypeNum)) %>%
    select(paramTypeNum, mean, se) -> summ2.noXcoefs
  summ2.pretty <- rbind(summ2.xcoefs, summ2.noXcoefs)
  
  #plot
  p <- ggplot(summ2.pretty, aes(x = paramTypeNum, 
                         y = mean, 
                         ymin = mean - se,
                         ymax = mean + se)) +
    geom_point() + 
    geom_errorbar(width=.1) + 
    coord_flip() +
    ylab("Converged parameters (%)") +
    xlab("Parameter type") + theme_bw() + ylim(c(0,100))
  p
  
  fileName <- paste0("output/boral_diagnostics/summarizeGeweke_byParamType_", modelID, "_", modelID2, ".pdf")
  ggsave(file = fileName, width = 3.5, height = 3)
  
}


#---------------------------------------------------------#
# plot mcmc chains

makeMCMCframe_allruns <- function(fit.list, test){
  
  #runType (used to make the iteration sequence)
  if(test == T){
    runType <- "test"
  }else{
    runType <- "real"
  }
  
  #extract mcmc obj and turn it into a df
  mcmc.list <- lapply(fit.list, function(x){
    mcmc.obj <- x$mcmc.obj
    mcmc.df <- data.frame(mcmc.obj[[1]])
    return(mcmc.df)
  })
  mcmc.df <- list_to_df(mcmc.list)
  
  # create sample iteration sequence
  mcmc_controls <- set_mcmc_controls(seed = 1, runType = runType)
  chain <- 1:mcmc_controls$n.iteration
  thin <- mcmc_controls$n.thin
  n.samples <- length(chain)/thin #sample every nth iteration in the chain
  samp <- list()
  for(i in 1:n.samples){
    samp[[i]] <- thin * i
  }
  samp <- unlist(samp)
  samp.afterburnin<-samp[samp > mcmc_controls$n.burnin]
  
  # extent sequence by the number of runs and add to the mcmc.df
  numRuns <- length(unique(mcmc.df$source))
  mcmc.df$iteration <- rep(samp.afterburnin, numRuns)
  
  # gather param cols
  mcmc.df %>%
    gather(key = "param", value = "estimate", -c(iteration, source)) -> mcmc.df.l
  
  # remove the deviance
  mcmc.df.l %>%
    filter(param != "deviance") -> mcmc.df.l
  
  return(mcmc.df.l)
  
}

extractMCMCfails <- function(mcmc.df.l, geweke.df.ann){
  
  #make the format of geweke.df.ann$param match mcmc.df.l$param by replacing all special characters with "."
  param.vec <- gsub("[", ".", geweke.df.ann$param, fixed=T)
  param.vec <- gsub(",", ".", param.vec, fixed=T)
  param.vec <- gsub("]", ".", param.vec, fixed=T)
  geweke.df.ann$param <- param.vec
  sum(unique(mcmc.df.l$param) != unique(geweke.df.ann$param)) # if 0, then these all match and that's awesome
  
  #join the mcmc and geweke dataframes
  mcmc.df.l %>%
    left_join(geweke.df.ann) -> mcmc.df.l.ann
  
  #extract mcmc fails
  mcmc.df.l.ann %>%
    filter(converg == "no") -> mcmc.fails
  
  return(mcmc.fails)
  
}

plotMCMC <- function(mcmc.df.l, modelLVonly){
  
  if(modelLVonly == TRUE){
    modelID <- "LV-only"
  }else{
    modelID <- "Traits-and-LVs"
  }
  
  mytheme <- make_ggplot_theme()
  
  #############################
  #lv.coefs
  mcmc.df.l %>%
    filter(grepl("lv.coefs", param)) %>%
    separate(param, into=c("drop1","drop2","OTUnum","paramNum", "drop3")) %>%
    select(source, iteration, OTUnum, paramNum, estimate) %>%
    mutate(param = "lv.coefs") %>%
    mutate(paramTypeNum = paste0(param, paramNum))-> lv.coefs.df
  
  p.lvcoefs <- ggplot(lv.coefs.df, aes(x = iteration, y = estimate, color = OTUnum)) + 
    geom_line() + 
    facet_grid(source ~ paramTypeNum, scales = "free") + 
    guides(color = F) +
    ggtitle("LV coef parameters, colors = OTUs") + 
    mytheme + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  #############################
  #lvs
  mcmc.df.l %>%
    filter(grepl("lvs", param)) %>%
    separate(param, into=c("drop1","SampNum","paramNum","drop2")) %>%
    select(source, iteration, SampNum, paramNum, estimate) %>%
    mutate(param = "lvs") %>%
    mutate(paramTypeNum = paste0(param, paramNum))-> lvs.df
  
  p.lvs <- ggplot(lvs.df, aes(x = iteration, y = estimate, color = SampNum)) + 
    geom_line() + 
    facet_grid(source ~ paramTypeNum, scales = "free") + 
    guides(color = F) +
    ggtitle("LVs, colors = samples") + 
    mytheme + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  #############################
  #row coefs
  #row coefs ID 1
  mcmc.df.l %>%
    filter(grepl("row.coefs.ID1", param)) %>%
    separate(param, into=c("drop1","drop2","drop3","paramNum","drop4")) %>%
    select(source, iteration, paramNum, estimate) %>%
    mutate(param = "rowID1") %>%
    mutate(OTUnum = NA) -> rowid1.df
  #row coefs ID 2
  mcmc.df.l %>%
    filter(grepl("row.coefs.ID2", param)) %>%
    separate(param, into=c("drop1","drop2","drop3","paramNum","drop4")) %>%
    select(source, iteration, paramNum, estimate) %>%
    mutate(param = "rowID2") %>%
    mutate(OTUnum = NA) -> rowid2.df
  #all row-related params
  rowid1.df %>%
    rbind(rowid2.df) -> tmp
  p.rows <- ggplot(tmp, aes(x = iteration, y = estimate, color = paramNum)) + 
    geom_line() + 
    facet_grid(source ~ param, scales = "free") + 
    guides(color = F) +
    ggtitle("Row coefs, colors = samples, sites") + 
    mytheme + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  #save all the plots
  l <- list(p.lvcoefs, p.lvs, p.rows)
  
  if(modelID == "Traits-and-LVs"){
    #############################
    #X coefs
    mcmc.df.l %>%
      filter(grepl("X.coefs", param)) %>%
      separate(param, into=c("drop1", "drop2","OTUnum","paramNum","drop3")) %>%
      select(source, iteration, OTUnum, paramNum, estimate) %>%
      mutate(param = "Xcoefs") -> Xcoefs.df
    XLEVELS <- unique(Xcoefs.df$paramNum)
    # first 6 X coefs
    Xcoefs.df %>%
      filter(paramNum %in% XLEVELS[1:6]) -> plot.df
    p.xcoefsA <- ggplot(plot.df, aes(x = iteration, y = estimate, color = OTUnum)) + 
      geom_line() + 
      facet_grid(source ~ paramNum, scales = "free") + 
      guides(color = F) +
      ggtitle("Xcoef parameters") + 
      mytheme + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    # last 5 X coefs
    Xcoefs.df %>%
      filter(paramNum %in% XLEVELS[7:11]) -> plot.df
    p.xcoefsB <- ggplot(plot.df, aes(x = iteration, y = estimate, color = OTUnum)) + 
      geom_line() + 
      facet_grid(source ~ paramNum, scales = "free") + 
      guides(color = F) +
      ggtitle("Xcoef parameters") + 
      mytheme + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    
    #save all the plots
    l <- list(p.lvcoefs, p.lvs, p.rows, p.xcoefsA, p.xcoefsB)
  }
  
  #save all the plots
  ml <- marrangeGrob(l, nrow = 1, ncol = 1)
  fileName <- paste0("output/boral_diagnostics/plotMCMC_", modelID, ".pdf")
  ggsave(file = fileName, ml)
  
}

makeFailplot_OTU <- function(mcmc.fails, summ, title, curr.paramType, curr.paramNum, selection){
  
  #make plot title
  summ %>%
    filter(paramTypeNum == paste0(curr.paramType, curr.paramNum)) -> curr.num
  plottitle <- paste0(title, ", ", curr.num$failedOTUs, " failed parameters")
  
  #make plot df
  mcmc.fails %>%
      filter(paramType == curr.paramType) %>%
      filter(paramNum2 == curr.paramNum) -> plot.df
  
  #subset failed OTUs
  plot.df %>%
    filter(OTUId_ann %in% unique(plot.df$OTUId_ann)[selection]) -> plot.df
  
  #plot
  mytheme <- make_ggplot_theme()
  p <- ggplot(plot.df, aes(x=iteration, y=estimate, color=source)) +
    geom_line() +
    facet_wrap(~OTUId_ann) + 
    ggtitle(plottitle) +
    mytheme + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_color_discrete(name = "Model run")
  
}

makeFailplot_samp <- function(mcmc.fails, summ, title, curr.paramType, curr.paramNum, row, selection){
  
  #make plot title
  summ %>%
    filter(paramTypeNum == paste0(curr.paramType, curr.paramNum)) -> curr.num
  plottitle <- paste0(title, ", ", curr.num$failedSamps, " failed parameters")
  
  #make plot df
  if(row == FALSE){
    mcmc.fails %>%
      filter(paramType == curr.paramType) %>%
      filter(paramNum2 == curr.paramNum) -> plot.df
  }else{
    mcmc.fails %>%
      filter(paramType == curr.paramType) -> plot.df
  }
  
  plot.df %>%
    filter(seq_sampName %in% unique(plot.df$seq_sampName)[selection]) -> plot.df
  
  #plot
  mytheme <- make_ggplot_theme()
  p <- ggplot(plot.df, aes(x=iteration, y=estimate, color=source)) +
    geom_line() +
    facet_wrap(~seq_sampName) + 
    ggtitle(plottitle) +
    mytheme + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_color_discrete(name = "Model run")
  
  return(p)
  
}

makeFailplot_rowID2 <- function(mcmc.fails, title, curr.paramType){
  
  #make plot title
  mcmc.fails %>%
    filter(paramType == curr.paramType) %>%
    summarize(failedParams = length(unique(paramNum1))) -> curr.num
  plottitle <- paste0(title, ", ", curr.num$failedParams, " failed parameters")
  
  #make plot df
  mcmc.fails %>%
    filter(paramType == curr.paramType) -> plot.df
  
  #plot
  mytheme <- make_ggplot_theme()
  p <- ggplot(plot.df, aes(x=iteration, y=estimate, color=source)) +
    geom_line() +
    facet_wrap(~paramNum1) + 
    ggtitle(plottitle) +
    mytheme + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_color_discrete(name = "Model run")
  
  return(p)
  
}

plotMCMC_fails <- function(mcmc.fails, modelLVonly){
  
  if(modelLVonly == TRUE){
    modelID <- "LV-only"
  }else{
    modelID <- "Traits-and-LVs"
  }
  
  mcmc.fails %>%
    group_by(paramTypeNum) %>%
    summarize(failedOTUs = length(unique(OTUId)),
              failedSamps = length(unique(sampNum))) -> summ
  
  selection <- 1:25
  
  # BY OTU
  #lv.coefs1
  p.lvcoef1 <- makeFailplot_OTU(mcmc.fails = mcmc.fails, summ = summ, selection = selection, 
                            title = "lv.coefs1 (OTU-specific intercepts)", 
                            curr.paramType = "lv.coefs", curr.paramNum = "1")
  
  #lv.coefs2
  p.lvcoef2 <- makeFailplot_OTU(mcmc.fails = mcmc.fails, summ = summ, selection = selection,
                            title = "lv.coefs2 (LV1 coefficients)", 
                            curr.paramType = "lv.coefs", curr.paramNum = "2")
  
  #lv.coefs3
  p.lvcoef3 <- makeFailplot_OTU(mcmc.fails = mcmc.fails, summ = summ, selection = selection,
                            title = "lv.coefs3 (LV2 coefficients)", 
                            curr.paramType = "lv.coefs", curr.paramNum = "3")
  
  #lv.coefs4
  p.lvcoef4 <- makeFailplot_OTU(mcmc.fails = mcmc.fails, summ = summ, selection = selection,
                            title = "lv.coefs4 (Dispersion)", 
                            curr.paramType = "lv.coefs", curr.paramNum = "4")
  
  # BY SAMP
  #lv1
  p.lv1 <- makeFailplot_samp(mcmc.fails = mcmc.fails, summ = summ, selection = selection,
                             title = "LV 1",
                             curr.paramType = "lvs", curr.paramNum = "1",
                             row = FALSE)
  
  #lv2
  p.lv2 <- makeFailplot_samp(mcmc.fails = mcmc.fails, summ = summ, selection = selection,
                             title = "LV 2",
                             curr.paramType = "lvs", curr.paramNum = "2",
                             row = FALSE)
  
  #row.coefs.ID1
  p.rowID1 <- makeFailplot_samp(mcmc.fails = mcmc.fails, summ = summ, selection = selection,
                                title = "RowID 1 (sample)", 
                                curr.paramType = "row.coefs.ID1", curr.paramNum = "NA",
                                row = TRUE)
  
  #BY SITE
  #row.coefs.ID2
  p.rowID2 <- makeFailplot_rowID2(mcmc.fails = mcmc.fails, 
                                  title = "RowID 2 (site)",
                                  curr.paramType = "row.coefs.ID2")
  
  #save all the plots
  l <- list(p.lvcoef1, p.lvcoef2, p.lvcoef3, p.lvcoef4, 
            p.lv1, p.lv2,
            p.rowID1, p.rowID2)
  ml <- marrangeGrob(l, nrow = 1, ncol = 1)
  fileName <- paste0("output/boral_diagnostics/plotMCMC_fails_", modelID, ".pdf")
  ggsave(file = fileName, ml)
  
}


#---------------------------------------------------------#
# How different are the mean parameter estimates that do and don't converge across runs?

little.plotfxn <- function(tmp){
  
  p <- ggplot(tmp, aes(x = reorder(OTUId_ann, desc(diff)), 
                       y = X, 
                       color=converg, shape = source)) +
    geom_point() +
    geom_errorbar(aes(ymin=X.lower, ymax=X.upper)) +
    coord_flip() +
    facet_grid(~term, scales = "free") +
    geom_hline(yintercept = 0, linetype = 2) +
    xlab("OTU (ordered by estimate difference)") + ylab("Estimate")
  
  return(p)
  
}

compareXestim_converg <- function(geweke.df.ann, fit.list, complete_subset_list, allXs){
  
  if(allXs == TRUE){
    modelID2 <- "allX"
    boral_dataobjs <- make_boral_dataobjs_allX(complete_subset_list)
  }else{
    modelID2 <- "selectX"
    boral_dataobjs <- make_boral_dataobjs_selectX(complete_subset_list)
  }
  
  #link the order of the covariates with their names
  X <- boral_dataobjs[["X"]]
  x.indx <- data.frame(term = colnames(X),
                       paramNum2 = as.character(seq(1:length(colnames(X)))))
  
  #extract X coef mean values
  xmean.list <- lapply(fit.list, function(x) x$Xcoefs.df)
  xmean.df <- list_to_df(xmean.list)
  xmean.df %>%
    left_join(x.indx) -> xmean.df
  
  #calc geweke.df.ann, filter for Xcoefs and merge with Xcoefs
  geweke.df.ann %>%
    filter(paramType == "X.coefs") %>%
    select(source, paramTypeNum, paramNum2, OTUId, OTUId_ann, numSamps_present, converg) %>%
    left_join(xmean.df) -> x.converg.df
  
  #for each X coef...
  XCOEF <- unique(x.converg.df$paramTypeNum)
  plotdf.list <- list()
  for(i in 1:length(XCOEF)){
    
    # filter by XCOEF
    x.converg.df %>%
      filter(paramTypeNum == XCOEF[i]) -> curr.x.converg.df
    
    #filter for the parameters that have runs where they converged/not
    curr.x.converg.df %>%
      group_by(OTUId_ann) %>%
      summarize(numCells = length(converg),
                numConv = length(unique(converg))) %>%
      filter(numConv == 2) -> curr.params
    curr.x.converg.df %>%
      filter(OTUId_ann %in% curr.params$OTUId_ann) -> plot.df
    
    #find the mean difference between the converged and not converged estimate
    plot.df %>%
      group_by(converg, OTUId_ann) %>%
      summarize(mean.X = mean(X)) %>%
      spread(key = converg, value = mean.X) %>%
      rename('mean.X.no'='no',
             'mean.X.yes' = 'yes') %>%
      mutate(diff = abs(mean.X.no - mean.X.yes)) %>%
      select(OTUId_ann, diff) -> diff.indx
    plot.df %>%
      left_join(diff.indx) -> plot.df
    
    plotdf.list[[i]]<-plot.df
    
  }
  l <- lapply(plotdf.list, little.plotfxn)
  #l[[1]] # density coefs
  
  #print plots
  # split the terms across separate pages
  ml <- marrangeGrob(l, nrow = 1, ncol = 1)
  fileName <- paste0("output/boral_diagnostics/compareXestim_converg_", modelID2, ".pdf")
  ggsave(file = fileName, ml)
  
}


#---------------------------------------------------------#
# residuals

plot_boral_residuals <- function(fit.list, modelLVonly, allXs){
  
  if(modelLVonly == TRUE){
    modelID <- "LV-only"
  }else{
    modelID <- "Traits-and-LVs"
  }
  
  if(allXs == TRUE){
    modelID2 <- "allX"
  }else{
    modelID2 <- "selectX"
  }
  
  resid.data <- lapply(fit.list, function(x) x$resid.data)
  resid.df <- list_to_df(resid.data)
  
  #plot
  mytheme <- make_ggplot_theme()
  p <- ggplot(data = resid.df, aes(x = log(fitted.vals), y = residuals, color = OTUId)) +
    geom_point(size = 0.5) + 
    geom_hline(yintercept = 0, linetype = 2) +
    facet_wrap(~source) +
    xlab("Linear predictor") + ylab("Dunn-Smyth residuals") +
    guides(color = FALSE) + mytheme
  p
  fileName <- paste0("output/boral_diagnostics/residuals_", modelID, "_", modelID2, ".pdf")
  ggsave(file = fileName, width = 6, height = 3)
  
  #plot(fit) #default residual plots from the fitted obj
}


#---------------------------------------------------------#
# ordination plots and the trace of the residual covariance matrix

plot_boral_ordinations <- function(fit.list, fit.lvOnly.list, allXs){
  
  if(allXs == TRUE){
    modelID2 <- "allX"
  }else{
    modelID2 <- "selectX"
  }
  
  #extract scaled lv coordinates
  scaled.lvs.lvonly <- lapply(fit.lvOnly.list, function(x) x$lvsplot.data$scaled.lvs)
  scaled.lvs.df.lvonly <- list_to_df(scaled.lvs.lvonly)
  scaled.lvs.df.lvonly$modelID <- "LV-only"
  scaled.lvs <- lapply(fit.list, function(x) x$lvsplot.data$scaled.lvs)
  scaled.lvs.df <- list_to_df(scaled.lvs)
  scaled.lvs.df$modelID <- "Traits-and-LVs"
  scaled.lvs.all <- full_join(scaled.lvs.df.lvonly, scaled.lvs.df)
  
  #plot and save
  mytheme <- make_ggplot_theme()
  p <- ggplot(scaled.lvs.all, aes(x = LV1, y = LV2)) + 
    geom_point(size=1) + coord_fixed() +
    facet_grid(modelID~source) + mytheme
  p
  fileName <- paste0("output/boral_diagnostics/ordinations_",modelID2, ".pdf")
  ggsave(file=fileName, height=6, width=12)
  
}

percExplained_byXvars <- function(fit.list, fit.lvOnly.list){
  
  trace.lvonly <- unlist(lapply(fit.lvOnly.list, function(x) x$rescor.trace))
  trace <- unlist(lapply(fit.list, function(x) x$rescor.trace))
  
  envcor <- ((mean(trace.lvonly) - mean(trace)) / mean(trace.lvonly)) *100
  # enviromental covariates can account for X% of the covariation between OTUs
  
  return(envcor)
  
}



