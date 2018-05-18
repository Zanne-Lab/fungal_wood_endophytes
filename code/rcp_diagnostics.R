# RCP diagnostics
#######


plot_numOfgroups_diagnostics <- function(rcp_output_results){
  
  # minnpp is the minimum value when taking the sum of membership probabilities for each cluster (very small values = empty groups = bad)
  # red line in lower right panel is a guide for even numbered groups
  
  rcp_output_results %>%
    mutate(minpp.goodfit = ifelse(minpp > 1, TRUE, FALSE)) -> plot.df 
  df <- data.frame(yy = (101/1:30), xx = 1:30) # yy is number of samples/ 1:max number of groups
  xlabname <- "Number of groups"
  
  mytheme <- make_ggplot_theme()
  p.bic <- ggplot(plot.df, aes(x = nRCP, y = BIC)) + 
    geom_point(pch = 16, alpha = .8) + geom_vline(aes(xintercept = 5), linetype=2) +
    xlab(xlabname) + mytheme
  p.aic <- ggplot(plot.df, aes(x = nRCP, y = AIC)) +
    geom_point(pch = 16, alpha = .8) + geom_vline(aes(xintercept = 5), linetype=2) +
    xlab(xlabname) + mytheme
  p.loglik <- ggplot(plot.df, aes(x = nRCP, y = logLik)) +
    geom_point(pch = 16, alpha = .8) + geom_vline(aes(xintercept = 5), linetype=2) +
    ylab("Penalized log-liklihood") + xlab(xlabname) + mytheme
  p.minpp <- ggplot(plot.df, aes(x = nRCP, y = minpp, shape = minpp.goodfit)) +
    geom_point(alpha = .8) + geom_vline(aes(xintercept = 5), linetype=2) +
    ylab("Sum of posterior membership probability") + xlab(xlabname) + mytheme +
    scale_shape_manual(name = "", values = c(1, 16), 
                       labels = c("Empty clusters","Good fit")) +
    theme(legend.position = c(0.7, 0.7))
  
  pdf(file = "output/rcp_diagnostics/numOfgroups_diagnostics.pdf", width = 6, height = 6)
  grid.arrange(p.bic + ggtitle("a"), 
               p.aic + ggtitle("b"), 
               p.loglik + ggtitle("c"), 
               p.minpp + ggtitle("d"), 
               ncol = 2)
  dev.off()
  
}

plot_rcp_residuals <- function(fm_rcp){
  pdf(file = "output/rcp_diagnostics/residuals.pdf", width = 6, height = 3)
  plot(fm_rcp)
  dev.off()
}

plot_rcpVSenv <- function(complete_subset_list, fm_rcp){
  
  covariates <- data.frame(complete_subset_list[["covariates"]])
  post_probs <- data.frame(fm_rcp$postProbs)
  names(post_probs) <- fm_rcp$names$RCPs
  covar_names <- c("Mn", "Fe", "Zn", "C", "P", "density")
  
  pdf(file = "output/rcp_diagnostics/rcp_preds-vs-environment.pdf", height = 6, width = 5)
  par(mfcol=c(3,2))
  for (rcp in names(post_probs)) {
    # plot(1, type="n", axes=F, xlab="", ylab="")
    # text(1,1,labels = rcp)
    for (covar in covar_names) {
      if (covar == "Mn") {
        plot(post_probs[,rcp] ~ covariates[,covar], 
             ylab = "", xlab = covar, main = rcp,
             col = alpha("black", 0.1), pch = 16)
      } else {
        plot(post_probs[,rcp] ~ covariates[,covar], 
             ylab = "", xlab = covar,
             col = alpha("black", 0.1), pch = 16)
      }
    }
  }
  dev.off()
  
}

HardClust <- function(postProbs){
  postProbs = data.frame(round(postProbs, 10))
  postProbs.rowMax = apply(postProbs,1,max)
  RCPclusters = integer(nrow(postProbs)) # pre-allocate vector
  for (i in 1:nrow(postProbs)){
    RCPclusters[i] = which(postProbs[i,]==postProbs.rowMax[i])
  }
  return(RCPclusters)
}

compare_rcp_to_boral <- function(fit.list, fit.lvOnly.list, fm_rcp){
  
  #just use the first run from each of the boral models
  fit.lvonly<- fit.lvOnly.list[[1]]
  fit.lvenv<- fit.list[[1]]
  
  #hard-cluster the samples
  rcp <- factor(HardClust(fm_rcp$postProbs)) ## function that returns hard cluster membership
  
  #make plotting dataframe
  rcp_boral <- data.frame(lv_only1 = fit.lvonly$lvmean.data[,1],
                          lv_only2 = fit.lvonly$lvmean.data[,2],
                          lv_env1 = fit.lvenv$lvmean.data[,1],
                          lv_env2 = fit.lvenv$lvmean.data[,2],
                          rcp = rcp)
  
  # plot 
  p1 <- ggplot(data = rcp_boral, aes(x=lv_only1, y=lv_only2, color=rcp)) + 
    geom_point(size=3) + 
    scale_colour_brewer("RCP \ncluster", palette = "Paired") +
    scale_x_continuous(limits = c(-3,1)) + # adjust as needed
    ggtitle("Boral ordination - latent variable only") +
    theme_bw()
  
  p2 <- ggplot(data = rcp_boral, aes(x=lv_env1, y=lv_env2, colour=rcp)) + 
    geom_point(size=3) + 
    scale_colour_brewer("RCP \ncluster", palette = "Paired") +
    ggtitle("Boral ordination - latent variable + env. covariates") +
    theme_bw()
  
  pdf(file = "output/rcp_diagnostics/plot_boralVSrcp.pdf", height = 6, width = 8)
  grid.arrange(p1, p2, ncol = 2)
  dev.off()
  
}
