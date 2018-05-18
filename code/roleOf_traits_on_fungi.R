# investigate the role of wood traits on fungal OTUs with vegan::capscale() and mvabund::manyglm()
#######

#---------------------------------------------------------#
# trait correlations

make_trait_cor_plot <- function(complete_subset_list){
  
  require(GGally)
  
  covariates <- complete_subset_list[['covariates']]
  covariates %>%
    select(C, waterperc, density, barkthick) -> plot.df
  
  ggpairs(plot.df)
  ggsave(file = "output/roleOfTraits/trait_cor.pdf", width = 6, height = 6)
  
}

#---------------------------------------------------------#
# distance-based methods

ggbiplot_custom <- function(mod.obj, envfit.obj, seqSamples){

  # collect site scores from mod.obj
  scrs <- as.data.frame(scores(mod.obj, display = "sites"))
  scrs <- cbind(scrs, seq_sampName = row.names(scrs))
  scrs %>%
    left_join(seqSamples) -> scrs
  
  # collect vector scores from envfit.obj
  spp.scrs <- as.data.frame(scores(envfit.obj, display = "vectors"))
  spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs), pval = envfit.obj$vectors$pvals)
  
  # spp.scrs %>%
  #   mutate(MDS1_nudgeRight = MDS1 + 0.05) %>%
  #   mutate(MDS1_nudgeLeft = MDS1 - 0.05) %>%
  #   mutate(MDS2_nudgeUp = MDS2 + 0.05) %>%
  #   mutate(MDS2_nudgeDown = MDS2 - 0.05) %>%
  #   mutate(MDS1_nudge = NA) %>%
  #   mutate(MDS2_nudge = NA) -> spp.scrs
  
  # upperRight <- c("N","barkthick","C")
  # upperLeft <- c("Zn","density")
  # lowerRight <- c("waterperc","P")
  # lowerLeft <- c("Fe","Mn","K","Ca")
  # 
  # select <- spp.scrs$Species %in% upperRight
  # spp.scrs[select, "MDS1_nudge"] <- spp.scrs[select, "MDS1_nudgeRight"]
  # spp.scrs[select, "MDS2_nudge"] <- spp.scrs[select, "MDS2_nudgeUp"]
  # 
  # select <- spp.scrs$Species %in% upperLeft
  # spp.scrs[select, "MDS1_nudge"] <- spp.scrs[select, "MDS1_nudgeLeft"]
  # spp.scrs[select, "MDS2_nudge"] <- spp.scrs[select, "MDS2_nudgeUp"]
  # 
  # select <- spp.scrs$Species %in% lowerRight
  # spp.scrs[select, "MDS1_nudge"] <- spp.scrs[select, "MDS1_nudgeRight"]
  # spp.scrs[select, "MDS2_nudge"] <- spp.scrs[select, "MDS2_nudgeDown"]
  # 
  # select <- spp.scrs$Species %in% lowerLeft
  # spp.scrs[select, "MDS1_nudge"] <- spp.scrs[select, "MDS1_nudgeLeft"]
  # spp.scrs[select, "MDS2_nudge"] <- spp.scrs[select, "MDS2_nudgeDown"]
  # 
  #spp.scrs %>%
  #  filter(pval < 0.05) -> signif.spp.scrs
  
  # plot
  mult <- 2 #multiplier for the arrows and text for envfit
  mult.text <- 2.3
  mytheme <- make_ggplot_theme()
  p <- ggplot(scrs) +
    geom_text(mapping = aes(x = MDS1, y = MDS2, label = code,
                             color = Family), size = 3) +
    coord_fixed() + ## need aspect ratio of 1!
    geom_segment(data = spp.scrs,
                 aes(x = 0, xend = mult*MDS1, y = 0, yend = mult*MDS2),
                 arrow = arrow(length = unit(0.2, "cm")), colour = 1) +
    geom_text(data = spp.scrs, 
              aes(x = mult.text*MDS1, y = mult.text*MDS2, 
                  label = Species), size = 4, fontface = "bold") +
    scale_color_viridis(name="Wood family", discrete=TRUE) +
    xlab("PCoA 1") + ylab("PCoA 2") + mytheme
  #    scale_shape_manual(name="Size class", values=c(16,1)) 
  
  return(p)
  
}

make_capscale_mods <- function(complete_subset_list){
  
  require(vegan)
  
  # identify the otu matrix and wood traits
  mat.otu <- complete_subset_list[["otus.trimmed"]]
  # traits agg by sp, stem size
  xVars <- complete_subset_list[["covariates"]]
  xVars %>%
    select(size, density, barkthick, waterperc, P, K, Ca, Mn, Fe, Zn, N, C) -> envVars
  
  # transform otu tab
  mat.otu.t<-decostand(mat.otu, 'hellinger')
  
  # unconstrained mod
  mod0 <- capscale(mat.otu.t ~ 1, data = envVars, distance='bray')
  
  # constrained mod
  modfull <- capscale(mat.otu.t ~ ., data=envVars, distance='bray')
  
  mod.list<-list(mod0 = mod0, modfull = modfull, envVars = envVars)
  return(mod.list)
  
}

plot_envfit_on_unconstrainedOrd <- function(complete_subset_list, seqSamples){
  
  # extract the unconstrained mod and covariates
  mod.list <- make_capscale_mods(complete_subset_list)
  mod0 <- mod.list[['mod0']]
  envVars <- mod.list[['envVars']]
  
  # do envfit
  mod0.envfit<-envfit(mod0, envVars)
  
  # make plot
  p <- ggbiplot_custom(mod.obj = mod0, envfit.obj = mod0.envfit, seqSamples = seqSamples)  
  p
  ggsave(file="output/roleOfTraits/unconstrained_envfit.pdf", width = 8, height = 5)
  
}

summary_bestmodel_constrainedOrd <- function(complete_subset_list){
  
  # extract the unconstrained and constrained mod
  mod.list <- make_capscale_mods(complete_subset_list)
  mod0 <- mod.list[['mod0']]
  modfull <- mod.list[['modfull']]
  
  # model selection
  modstep <- ordistep(mod0, scope=formula(modfull))
  
  # summary table for the best model
  df <- data.frame(term=row.names(modstep$anova), modstep$anova)
  df %>%
    separate(term, into=c("drop","term")) %>%
    select(-drop) %>%
    rename('pval'=`Pr..F.`) -> df
  
  write.csv(df, file = "output/roleOfTraits/constrained_bestmodelsummary.csv")
  
}


#---------------------------------------------------------#
# mvabund methods

fit_mvabund<-function(complete_subset_list, minPerc){
  
  require(mvabund)
  #The mvabund approach improves power across a range of species with different variances and includes an assumption of a mean-variance relationship. 
  #It does this by fitting a single generalised linear model (GLM) to each response variable with a common set of predictor variables
  
  # identify the otu matrix and wood traits, scale the continuous wood traits
  mat.otu <- complete_subset_list[["otus.trimmed"]]
  xVars <- complete_subset_list[["covariates"]]
  xVars %>%
    select(seq_sampName, species, site, size,
           density, barkthick, waterperc, P, K, Ca, Mn, Fe, Zn, N, C) -> envVars
  
  # specify data for mvabund
  fung <- list(abund = mat.otu, x = envVars)
  Dat <- mvabund(fung$abund, row.names=row.names(fung$abund)) # so R knows to treat coverDat as multivariate abundance
  size <- factor(fung$x$size)
  density <- fung$x$density
  barkthick <- fung$x$barkthick
  waterperc <- fung$x$waterperc
  P <- fung$x$P
  K <- fung$x$K
  Ca <- fung$x$Ca
  Mn <- fung$x$Mn
  Fe <- fung$x$Fe
  Zn <- fung$x$Zn
  N <- fung$x$N
  C <- fung$x$C
  
  # specify mvabund models
  # default family = "negative binomial", assumes a quadratic mean-variance relationship and a log-linear relationship between the response variables and any continuous variables
  ft.base <- manyglm(Dat ~ size)
  ft.m.density <- manyglm(Dat ~ size + density)
  ft.m.barkthick <- manyglm(Dat ~  size + barkthick)
  ft.m.waterperc <- manyglm(Dat ~ size + waterperc)
  ft.m.P <- manyglm(Dat ~ size +  P)
  ft.m.K <- manyglm(Dat ~ size +  K)
  ft.m.Ca <- manyglm(Dat ~ size + Ca)
  ft.m.Mn <- manyglm(Dat ~ size + Mn)
  ft.m.Fe <- manyglm(Dat ~ size + Fe)
  ft.m.Zn <- manyglm(Dat ~ size + Zn)
  ft.m.N <- manyglm(Dat ~ size + N)
  ft.m.C <- manyglm(Dat ~ size + C)
  mod.m.list <-list(base = ft.base, density = ft.m.density, barkthick = ft.m.barkthick, waterperc = ft.m.waterperc,
                    P = ft.m.P, K = ft.m.K, Ca = ft.m.Ca, Mn = ft.m.Mn, Fe = ft.m.Fe, Zn = ft.m.Zn, N = ft.m.N, C = ft.m.C, Dat = Dat)
  return(mod.m.list)
}

meanvariance_plot <- function(mod.m.list){
  
  Dat <- mod.m.list[['Dat']]
  
  pdf(file = "output/roleOfTraits/meanVariance.pdf", width = 5, height = 5)
  meanvar.plot(Dat, xlab="Mean OTU abundance", ylab="Variance in OTU abundance") # with log scale
  dev.off()
  
}

mvabund_residualvFitted_plot <- function(mod){
  
  require(ggplot2)
  
  #plot(mod, which=1, cex=0.5, caption = "", xlab = "")
  #why the heck does the default plot look different?  It looks like linear predictors less than -6 have been cut out
  
  #extract linear predictors
  seq_sampName <- row.names(mod$linear.predictor)
  linearpred.mat <- mod$linear.predictor
  linearpred.df <- data.frame(seq_sampName, linearpred.mat)
  
  #extract residuals
  residual.mat<-residuals(mod)
  residual.df<-data.frame(seq_sampName, residual.mat)
  colnames(residual.df) <- colnames(linearpred.df)
  
  #make plotting df
  linearpred.df %>%
    gather(key="OTUId", value="linear.predictor", -seq_sampName) -> linearpred.df
  residual.df %>%
    gather(key="OTUId", value="residuals", -seq_sampName) -> residual.df
  linearpred.df %>%
    left_join(residual.df) -> resid.plot.df
  
  #plot
  mytheme <- make_ggplot_theme()
  p <- ggplot(data = resid.plot.df, aes(x = linear.predictor, y = residuals, color = OTUId)) +
    geom_point(size=0.5) + 
    geom_hline(yintercept = 0, linetype = 2) +
    xlab("Linear predictor") + ylab("Dunn-Smyth residuals") +
    guides(color=FALSE) + mytheme
  
  return(p)
  
}

check_mvabundAssumpts <- function(mod.m.list){
  
  require(gridExtra)
  
  mod.list<-mod.m.list[names(mod.m.list) != "Dat"] #strip off the data piece to get just the list of mvabund models
  p.list <- lapply(mod.list, mvabund_residualvFitted_plot)
  
  #add plot titles
  p.list.titled<-list()
  for(i in 1:length(p.list)){
    p.list.titled[[i]] <- p.list[[i]] + ggtitle(names(mod.list)[i])
  }
  
  #layout
  n <- length(p.list.titled)
  nCol <- floor(sqrt(n))
  
  #plot
  pdf(file = "output/roleOfTraits/mvabund_residVfitted.pdf")
  do.call("grid.arrange", c(p.list.titled, ncol=nCol))
  dev.off()
  
  #note that different colors represent different OTUs
  
}

write_mvabund_aicTable<-function(mod.m.list){
  
  #strip off the data piece to get just the list of mvabund models
  mod.list<-mod.m.list[names(mod.m.list) != "Dat"] 
  
  #calculate the total AIC for each model
  aic.list <- lapply(mod.list, function(x){
    sum(AIC(x)) # sum AIC values for each OTUId
  })
  
  #calculate deltaAIC = base AIC - candidate AIC
  baseAIC <- aic.list[['base']]
  aic.df <- data.frame(modelName = names(unlist(aic.list)), 
             AIC = unlist(aic.list))
  aic.df %>%
    mutate(deltaAIC = baseAIC - AIC) %>%
    arrange(desc(deltaAIC)) -> aic.df
  
  write.csv(aic.df, file = "output/roleOfTraits/mvabundAICs.csv")
  
}
