# boral_diagnostics fxns that are used in testing code

makeMCMCframe_allruns_customMCMCcontrols <- function(fit.list, mcmc_controls){
  
  #extract mcmc obj and turn it into a df
  mcmc.list <- lapply(fit.list, function(x){
    mcmc.obj <- x$mcmc.obj
    mcmc.df <- data.frame(mcmc.obj[[1]])
    return(mcmc.df)
  })
  mcmc.df <- list_to_df(mcmc.list)
  
  # create sample iteration sequence
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

makeFailplot_OTU_test <- function(mcmc.fails, title, curr.paramType, curr.paramNum, selection){
  
  #make plot title
  plottitle <- paste0(title, ", failed parameters")
  
  #make plot df
  mcmc.fails %>%
    filter(paramType == curr.paramType) %>%
    filter(paramNum2 == curr.paramNum) -> plot.df
  
  #subset failed OTUs
  plot.df %>%
    filter(OTUId_ann %in% unique(plot.df$OTUId_ann)[selection]) -> plot.df
  
  #plot
  p <- ggplot(plot.df, aes(x=iteration, y=estimate, color=source)) +
    geom_line() +
    facet_grid(test~OTUId_ann) + 
    ggtitle(plottitle) +
    theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_color_discrete(name = "Model run")
  
  return(p)
  
}

makeFailplot_samp_test <- function(mcmc.fails, title, curr.paramType, curr.paramNum, row, selection){
  
  #make plot title
  plottitle <- paste0(title, ", failed parameters")
  
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
  p <- ggplot(plot.df, aes(x=iteration, y=estimate, color=source)) +
    geom_line() +
    facet_grid(test ~ seq_sampName) + 
    ggtitle(plottitle) +
    theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_color_discrete(name = "Model run")
  
  return(p)
  
}

makeFailplot_rowID2_test <- function(mcmc.fails, title, curr.paramType){
  
  #summarize
  mcmc.fails %>%
    filter(paramType == curr.paramType) %>%
    group_by(test) %>%
    summarize(failedParams = length(unique(paramNum1))) -> curr.num
  
  #make plot title
  plottitle <- paste0(title, ", failed parameters")
  
  #make plot df
  mcmc.fails %>%
    filter(paramType == curr.paramType) -> plot.df
  
  #plot
  p <- ggplot(plot.df, aes(x=iteration, y=estimate, color=source)) +
    geom_line() +
    facet_grid(test ~ paramNum1) + 
    ggtitle(plottitle) +
    theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_color_discrete(name = "Model run")
  
  return(p)
  
}

create_geweke_df_lowAnn <- function(fit.list){
  
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
  
  return(geweke.df.ann)
  
}

plot_cor_distributions_allruns_p.list <- function(fit.list, taxAndFunguild){
  
  cor.list <- lapply(fit.list, function(x) x$cor.df)
  cor.df <- list_to_df(cor.list)
  cor.df.ann <- annotate_cor_withOTUInfo(cor.df = cor.df, taxAndFunguild = taxAndFunguild)
  
  #plotting params
  mytheme <- make_ggplot_theme()
  colorVec <- cor_colors()
  names(colorVec) <- c("negative","positive","none")
  
  #signif factors
  cor.df.ann$enviro_signif <- recode(cor.df.ann$enviro_signif, no = "none")
  cor.df.ann$residual_signif <- recode(cor.df.ann$residual_signif, no = "none")
  signif.order <- c("none","positive","negative")
  cor.df.ann$enviro_signif <- factor(cor.df.ann$enviro_signif, levels = signif.order)
  cor.df.ann$residual_signif <- factor(cor.df.ann$residual_signif, levels = signif.order)
  
  #pair_haveGenusIDs
  cor.df.ann$pair_haveGenusIDs <- factor(cor.df.ann$pair_haveGenusIDs, levels = c("yes","no"))
  
  #enviro
  p.env<-ggplot(cor.df.ann, aes(x = enviro_cor, fill = enviro_signif, alpha = pair_haveGenusIDs)) + 
    geom_histogram() +
    geom_vline(xintercept = 0, linetype=2) +
    xlab("Shared environment correlation") + ylab("Frequency") +
    scale_fill_manual(name = "Signif.", values = colorVec) +
    scale_alpha_manual(name = "Genera IDs", values = c(0.5, 1)) +
    mytheme +
    facet_grid(~source) + guides(fill=FALSE, alpha=FALSE)
  #annotate("text", y=100, x=-0.5, label=paste(round(percNeg, digits=0),"%", sep="")) +
  #annotate("text", y=350, x=0.5, label=paste(round(percPos, digits=0),"%", sep="")) +
  
  #residual
  p.res<-ggplot(cor.df.ann, aes(x = residual_cor, fill = residual_signif, alpha = pair_haveGenusIDs)) + 
    geom_histogram() +
    geom_vline(xintercept = 0, linetype=2) +
    xlab("Residual correlation") + ylab("Frequency") +
    scale_fill_manual(name = "Signif.", values = colorVec) +
    scale_alpha_manual(name = "Genera IDs", values = c(0.5, 1)) +
    mytheme +
    facet_grid(~source) + guides(fill=FALSE, alpha=FALSE) +
    theme(axis.text.x = element_text(angle=90))
  #annotate("text", y=100, x=-0.5, label=paste(round(percNeg, digits=0),"%", sep="")) +
  #annotate("text", y=350, x=0.5, label=paste(round(percPos, digits=0),"%", sep="")) +
  
  p.list <- list(p.env = p.env, p.res = p.res)
  return(p.list)
  
}

barplot_cor_signif_p.list <- function(fit.list){
  
  cor.list <- lapply(fit.list, function(x) x$cor.df)
  cor.df <- list_to_df(cor.list)
  cor.df %>%
    mutate(otupair = paste0(otu1, otu2)) -> cor.df.ann
  
  # how many pairs in total?
  numOTUpairs <- length(unique(cor.df.ann$otupair))
  
  #summarize
  cor.df.ann %>%
    group_by(enviro_signif, source) %>%
    summarize(n_pairs = sum(!is.na(enviro_cor)),
              n_perc = (n_pairs/numOTUpairs) * 100 ,
              mean = round(mean(enviro_cor), digits = 2)) %>%
    mutate(corType = "enviro") %>%
    rename('sign'='enviro_signif') -> summ.enviro
  
  cor.df.ann %>%
    group_by(residual_signif, source) %>%
    summarize(n_pairs = sum(!is.na(residual_cor)),
              n_perc = (n_pairs/numOTUpairs) * 100 ,
              mean = round(mean(residual_cor), digits = 2)) %>%
    mutate(corType = "residual") %>%
    rename('sign'='residual_signif') -> summ.residual
  
  #plot
  plot.df <- rbind(summ.enviro, summ.residual)
  plot.df %>%
    filter(sign != "no") -> plot.df
  p.perc <- ggplot(plot.df, aes(x = sign, y = n_perc, label = mean)) +
    geom_point(pch = 1, size = 3) +
    facet_grid(~corType) +
    ylab("Percent of OTU pairs") + xlab("Correlation sign")
  
  return(p.perc)
  
}

# corType must be "enviro" or "residual"
cor.btwRuns_p.list <- function(fit.list, taxAndFunguild, corType){
  
  cor.list <- lapply(fit.list, function(x) x$cor.df)
  cor.df <- list_to_df(cor.list)
  cor.df.ann <- annotate_cor_withOTUInfo(cor.df = cor.df, taxAndFunguild = taxAndFunguild)
  cor.df.ann %>%
    mutate(otupair = paste0(otu1, otu2)) -> cor.df.ann
  
  # how many pairs in total?
  numOTUpairs <- length(unique(cor.df.ann$otupair))
  
  if(corType == "enviro"){
    cor.df.ann %>%
      select(source, otupair, enviro_cor, enviro_signif) %>%
      rename('cor'='enviro_cor',
             'signif'='enviro_signif') -> cor.df.ann
  }
  if(corType == "residual"){
    cor.df.ann %>%
      select(source, otupair, residual_cor, residual_signif) %>%
      rename('cor'='residual_cor',
             'signif'='residual_signif') -> cor.df.ann
    
  }
  
  # how many pairs change significance categories between runs?
  cor.df.ann %>%
    select(source, otupair, cor, signif) %>%
    group_by(otupair) %>%
    summarize(numSignifLevels = length(unique(signif)),
              signifLevels = paste(unique(signif), collapse = "_")) -> summ
  summ %>%
    filter(numSignifLevels !=1) -> tmp
  flipCat.otupairs <- unique(tmp$otupair)
  numFlip <- length(flipCat.otupairs)
  percFlip <- round((numFlip / numOTUpairs) * 100, digits = 0)
  
  # calculate the maximum difference between runs for otupairs and annotate the data
  cor.df.ann %>%
    select(source, otupair, cor) %>%
    group_by(otupair) %>%
    summarize(max = max(cor),
              min = min(cor)) %>%
    mutate(diff = abs(max - min)) %>%
    left_join(summ) %>%
    select(otupair, diff, numSignifLevels) -> diff.indx
  cor.df.ann %>%
    select(source, otupair, cor, signif) %>%
    left_join(diff.indx) -> plot.df
  
  # plot the distribution of differences for otupairs, color based on if the pair flips significance categories
  PAIR<-unique(plot.df$otupair)
  title <- paste(numFlip, " of ", numOTUpairs, 
                 " (", percFlip, "%) ",
                 "OTU pairs flip signif", sep ="")
  p.hist <- ggplot(plot.df, aes(x = diff, fill = factor(numSignifLevels))) +
    geom_histogram() + xlab("Max diff in correlation values across runs") +
    ggtitle(title) +
    guides(fill = F)
  #p.hist
  
  # plot the correlation values for OTU pairs where the difference between runs is greater that 0.5
  plot.df %>%
    filter(numSignifLevels != 1) %>%
    filter(diff > .3) -> plot.df.sub
  p.badones <- ggplot(plot.df.sub, aes(y = reorder(otupair, diff), x = cor)) +
    geom_point(aes(color = signif)) +
    geom_line() +
    geom_vline(aes(xintercept = 0), linetype = 2) +
    guides(color = F) +
    xlab("Correlation") + ylab("OTU pair")
  #p.badones
  
  p.list.list <- list(hist = p.hist, badones = p.badones)
  
  return(p.list.list)
  
}

