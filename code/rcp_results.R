# RCP results
#######

#---------------------------------------------------------#
# wood species confusion matrix

create_confusion <- function(postProbs, compare_to_data, compare_to_column) {
  
  # build binary matrix for categorical class being compared to
  compare_to <- model.matrix(as.formula(paste0("~0+",compare_to_column)), data=compare_to_data)
  # build matrix of expected shared sites
  sharedSites = t(postProbs) %*% compare_to
  sharedSites = round(sharedSites, 1)
  data.frame(sharedSites)
  
}

make_confus_table <- function(fm_rcp, complete_subset_list, seqSamples){
  
  covariates <- complete_subset_list[["covariates"]]
  # gives expected number of shared sites (accoutning for probabilistic membership) between RCPs and a another clustering variable
  # e.g. we expect 6 samples to be classified as RCP3 and come from the wood species alli and no samples that come from the wood species alli are expected to be classified as any other RCP group.
  rcp_species_confmat <- create_confusion(postProbs = fm_rcp$postProbs, 
                                          compare_to_data = covariates, 
                                          compare_to_column = "species")
  
  # make table pretty
  confmat.df <- data.frame(t(rcp_species_confmat))
  colnames(confmat.df) <- paste0("RCP", 1:5)
  species <- list_to_df(strsplit(row.names(confmat.df), "species"))[,2]
  confmat.df <- data.frame(species, confmat.df)
  sp.indx <- unique(seqSamples[,c("species","Binomial")])
  confmat.df %>%
    left_join(sp.indx) %>%
    select(Binomial, RCP1, RCP2, RCP3, RCP4, RCP5) %>%
    rename('Wood species' = Binomial) -> confmat.df
  
  return(confmat.df)
}

write_confus_table <- function(fm_rcp, complete_subset_list, seqSamples){
  confmat.df <- make_confus_table(fm_rcp = fm_rcp, complete_subset_list = complete_subset_list, seqSamples = seqSamples)
  write.csv(confmat.df, file="output/rcp_results/species_confmat.csv", row.names = FALSE)
}

plot_confus <- function(fm_rcp, complete_subset_list, seqSamples, zanneTree){
  
  # make confusion matrix
  confmat.df <- make_confus_table(fm_rcp = fm_rcp, complete_subset_list = complete_subset_list, seqSamples = seqSamples)
  colnames(confmat.df)[1] <- "species"
  # gives expected number of shared sites (accoutning for probabilistic membership) between RCPs and a another clustering variable
  
  # reformat df for ggtree
  row.names(confmat.df)<-sub(" ","_", confmat.df$species)
  confmat.df <- confmat.df[,-1] 
  graphing_df <- apply(confmat.df, 2,as.numeric)
  row.names(graphing_df) <- row.names(confmat.df)
  graphing_df[graphing_df==0] <- NA
  
  # color
  colorVec <- blueGradient_colors()
  
  # plot
  ggt <- ggtree(zanneTree) + geom_tiplab(size=3)
  gheatmap(ggt, graphing_df, 
           offset = 125, width=0.9, 
           colnames_position="top",
           font.size=2.5) + 
    scale_fill_gradient(na.value='white', low = colorVec["low"], high = colorVec["posBlue"])
  
  ggsave(file = "output/rcp_results/species_confmat.pdf", width = 6, height = 6)
}

plot_studyDesign <- function(seqSamples, zanneTree){
  
  #study design
  seqSamples %>%
    group_by(Binomial, size) %>%
    summarize(n = length(seq_sampName)) %>%
    spread(key = size, value = n) -> studydesign_df
  
  #pretty graphing df
  graphing_df<-data.frame(Small=studydesign_df$small, Large=studydesign_df$large)
  row.names(graphing_df)<-sub(" ","_",studydesign_df$Binomial)
  graphing_df[!is.na(graphing_df)] <- 1 # turn into presence/absence
  
  ggt <- ggtree(zanneTree) + geom_tiplab(size=4)
  gheatmap(ggt, graphing_df, 
                offset = 100, width=0.4, 
                colnames_position="top",
                font.size=4) +
    scale_fill_gradient(high = "black",low = "white", na.value = "white") +
    guides(fill = FALSE)
  
  ggsave(file = "output/rcp_results/species_studydesign.pdf", width = 6, height = 6)
  
}


#---------------------------------------------------------#
# partial effects plots

colMedian <- function(x) {
  meds <- numeric(ncol(x))
  for (i in 1:ncol(x)){
    meds[i] <- median(x[,i])
  }
  return(meds)
}

make_newdat <- function(covar.data, sizeLevel){
  
  # make a new data matrix
  new.dat <- matrix(NA, ncol=ncol(covar.data), nrow=1000)
  
  # get the col means for the continuous vars
  covar.data %>%
    select(-size) -> contin.covar.data
  contin.X.means <- colMeans(contin.covar.data)
  # get the col medians
  contin.X.meds <- colMedian(contin.covar.data)
  
  # add meds or means to the new data matrix?
  #lapply(X = contin.covar.data, FUN = hist)
  # means seems to make more sense here, medians would be odd
  
  #populate new.dat
  for (i in 1:length(contin.X.means)) {
    new.dat[,1+i] <- contin.X.means[i]
  }
  new.dat <- data.frame(new.dat)
  names(new.dat) <- names(covar.data)
  
  # make the polynomials - ensure names match original call
  # you could do this via a lapply call for the poly() calcs, and use model.covariate.vector with some tricky renaming,
  # but i like to see the varaibles laid out for some reason.
  new.dat <- data.frame(
    factor(sizeLevel, levels = c("small","large")), 
    predict(poly(covar.data$density, 2), newdata=new.dat$density),
    predict(poly(covar.data$barkthick, 2), newdata=new.dat$density),
    predict(poly(covar.data$waterperc, 2), newdata=new.dat$waterperc),
    predict(poly(covar.data$P, 2), newdata=new.dat$P),
    predict(poly(covar.data$K, 2), newdata=new.dat$K),
    predict(poly(covar.data$Ca, 2), newdata=new.dat$Ca),
    predict(poly(covar.data$Mn, 2), newdata=new.dat$Mn),
    predict(poly(covar.data$Fe, 2), newdata=new.dat$Fe),
    predict(poly(covar.data$Zn, 2), newdata=new.dat$Zn),
    predict(poly(covar.data$N, 2), newdata=new.dat$N),
    predict(poly(covar.data$C, 2), newdata=new.dat$C))
  
  names(new.dat) <- c(
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
  
  return(new.dat)
  
}

ParPlotMany.regimix <- function(fm, covar.data, newdata, variable, RCPs, 
                                xlab, ylab, legend, col, lty, title, group_names, xaxis, yaxis, ymtext1, ymtext2) {
    
    # identify the variable range (based on sizeLevel?)
    sizeLevel <- unique(newdata[,1])
    covar.data %>%
      filter(size == sizeLevel) -> covar.data.subset
    #var.range = seq(min(covar.data.subset[,variable]),max(covar.data.subset[,variable]),length.out=1000)
    var.range = seq(min(covar.data[,variable]),max(covar.data[,variable]),length.out=1000)
    
    # populate variable with values along the observed range instead of mean val
    newdata[,grep(variable, colnames(newdata))] = predict(poly(covar.data[,variable],2),newdata=var.range)
    
    # do prediction
    preds = predict.regimix(fm, newdata=newdata, nboot=0)
    
    if(legend == FALSE){
      
      # plot it
      matplot(x = var.range, 
              y = preds[,RCPs], 
              type = 'l', axes = FALSE,
              lwd = 3, 
              lty = lty, col = col,  
              ylab = "", xlab = "", ylim = c(0,1))
      
      title(main = title, adj = 0, line = 0.5, font.main = 2)
      rug(covar.data.subset[,variable]) # this part is important to subset by sizeLevel
      box()
      
      if(yaxis == TRUE){
        axis(side = 2) # yaxis
        mtext(text = ylab, side = 2, line = 3)
      }
      
      if(xaxis == TRUE){
        axis(side = 1) # xaxis
        mtext(text = xlab, side = 1, line = 3)
      }
      
      if(ymtext1 == TRUE){
        mtext(text = "Small stems", side = 2, font = 4, line = 5)
      }
      
      if(ymtext2 == TRUE){
        mtext(text = "Large stems", side = 2, font = 4, line = 5)
      }
      
    }else{
      
      # plot it
      matplot(x = var.range, 
              y = preds[,RCPs], 
              type = 'n', axes = FALSE,
              lwd = 3, 
              lty = lty, col = col,  
              ylab = "", xlab = "", ylim = c(0,1))
      
      legend(x = 30, y = 0.6, 
             legend = group_names, lty = lty,
             col = col, lwd = 3, horiz = FALSE)
      
    }
    
    
}

plot_partialEffects_select4 <- function(complete_subset_list, fm_rcp){
  
  # load covariate info
  rcp.objs <- make_rcp_dataobjs(complete_subset_list) # function written in 'rcp_modfit.R'
  covar.data <- rcp.objs$covar.data.forpredict
  
  # load data
  new.dat.small <- make_newdat(covar.data = covar.data, sizeLevel = "small")
  new.dat.large <- make_newdat(covar.data = covar.data, sizeLevel = "large")
  
  # plotting params
  
  # ... for RCPs
  col <- rcp_colors()
  lty <- rep(1, 5)
  group_names <- paste0("Group ",1:5)
  
  # ... for variables
  variable <- c("waterperc", "density", "C", "barkthick")
  xlab <- list("Water (%)", expression(Density ~ (g/cm^{3})), "Carbon (%)", "Bark Thickness (mm)")
  
  # ... for sizeLevel
  new.dat.list <- list(small = new.dat.small, large = new.dat.large)
  
  # ... for panels
  ylab <- "Membership probability"
  npanels <- length(new.dat.list) * length(variable)
  title.mat <- matrix(letters[1:npanels], nrow = length(new.dat.list), ncol = length(variable), byrow=T)
  legend.mat <- xaxis.mat <- yaxis.mat <- ymtext1.mat <- ymtext2.mat <- matrix(FALSE, nrow = length(new.dat.list), ncol = length(variable))
  xaxis.mat[2,] <- TRUE
  yaxis.mat[,1] <- TRUE
  ymtext1.mat[1,1] <- TRUE
  ymtext2.mat[2,1] <- TRUE
  
  # save
  pdf(file = "output/rcp_results/partialEffects_select4.pdf", width = 8.5, height = 5)
  par(mfcol=c(2,5), cex=.8) # tcl = -0.3
  par(oma=c(5,6,1,1)) # outer margins: bottom, left, top, right
  par(mar=c(0,0,1,0)+.5) # panel margins bottom, left, top, right
  
  for (i in 1:length(variable)) {
    
    for(k in 1:length(new.dat.list)){
      
      title <- title.mat[k,i]
      legend <- legend.mat[k,i]
      xaxis <- xaxis.mat[k,i]
      yaxis <- yaxis.mat[k,i]
      ymtext1 <- ymtext1.mat[k,i]
      ymtext2 <- ymtext2.mat[k,i]
      
      ParPlotMany.regimix(variable=variable[i], xlab=xlab[[i]],
                          fm=fm_rcp, 
                          newdata=new.dat.list[[k]],
                          covar.data=covar.data, 
                          RCPs=1:5, group_names=group_names, col=col, lty=lty, 
                          ylab=ylab, title=title, 
                          legend=legend, xaxis=xaxis, yaxis=yaxis, 
                          ymtext1=ymtext1, ymtext2=ymtext2)
      
    }
  }
  
  # add legend
  ParPlotMany.regimix(variable="waterperc", xlab=NULL,
                      fm=fm_rcp, 
                      newdata=new.dat.list[[1]],
                      covar.data=covar.data, 
                      RCPs=1:5, group_names=group_names, col=col, lty=lty, 
                      ylab=ylab, title=title, 
                      legend=TRUE, xaxis=xaxis, yaxis=yaxis, 
                      ymtext1=ymtext1, ymtext2=ymtext2)
  
  
  dev.off()
  
}

plot_partialEffect_theRest1 <- function(complete_subset_list, fm_rcp){
  
  # load covariate info
  rcp.objs <- make_rcp_dataobjs(complete_subset_list) # function written in 'rcp_modfit.R'
  covar.data <- rcp.objs$covar.data.forpredict
  
  # load data
  new.dat.small <- make_newdat(covar.data = covar.data, sizeLevel = "small")
  new.dat.large <- make_newdat(covar.data = covar.data, sizeLevel = "large")
  
  # plotting params
  
  # ... for RCPs
  col <- rcp_colors()
  lty <- rep(1, 5)
  group_names <- paste0("Group ",1:5)
  
  # ... for variables
  variable <- c("N","P","K")
  xlab <- list("Nitrogen (%)", "Phosphorus (%)", "Potassium (%)")
  
  # ... for sizeLevel
  new.dat.list <- list(small = new.dat.small, large = new.dat.large)
  
  # ... for panels
  ylab <- "Membership probability"
  npanels <- length(new.dat.list) * length(variable)
  title.mat <- matrix(letters[1:npanels], nrow = length(new.dat.list), ncol = length(variable), byrow=T)
  legend.mat <- xaxis.mat <- yaxis.mat <- ymtext1.mat <- ymtext2.mat <- matrix(FALSE, nrow = length(new.dat.list), ncol = length(variable))
  xaxis.mat[2,] <- TRUE
  yaxis.mat[,1] <- TRUE
  ymtext1.mat[1,1] <- TRUE
  ymtext2.mat[2,1] <- TRUE
  
  # save
  pdf(file = "output/rcp_results/partialEffects_theRest1.pdf", width = 7.5, height = 5)
  par(mfcol=c(2,4), cex=.8) # tcl = -0.3
  par(oma=c(5,6,1,1)) # outer margins: bottom, left, top, right
  par(mar=c(0,0,1,0)+.5) # panel margins bottom, left, top, right
  
  for (i in 1:length(variable)) {
    
    for(k in 1:length(new.dat.list)){
      
      title <- title.mat[k,i]
      legend <- legend.mat[k,i]
      xaxis <- xaxis.mat[k,i]
      yaxis <- yaxis.mat[k,i]
      ymtext1 <- ymtext1.mat[k,i]
      ymtext2 <- ymtext2.mat[k,i]
      
      ParPlotMany.regimix(variable=variable[i], xlab=xlab[[i]],
                          fm=fm_rcp, 
                          newdata=new.dat.list[[k]],
                          covar.data=covar.data, 
                          RCPs=1:5, group_names=group_names, col=col, lty=lty, 
                          ylab=ylab, title=title, 
                          legend=legend, xaxis=xaxis, yaxis=yaxis, 
                          ymtext1=ymtext1, ymtext2=ymtext2)
      
    }
  }
  
  # add legend
  ParPlotMany.regimix(variable="waterperc", xlab=NULL,
                      fm=fm_rcp, 
                      newdata=new.dat.list[[1]],
                      covar.data=covar.data, 
                      RCPs=1:5, group_names=group_names, col=col, lty=lty, 
                      ylab=ylab, title=title, 
                      legend=TRUE, xaxis=xaxis, yaxis=yaxis, 
                      ymtext1=ymtext1, ymtext2=ymtext2)
  
  dev.off()
  
}

plot_partialEffect_theRest2 <- function(complete_subset_list, fm_rcp){
  
  # load covariate info
  rcp.objs <- make_rcp_dataobjs(complete_subset_list) # function written in 'rcp_modfit.R'
  covar.data <- rcp.objs$covar.data.forpredict
  
  # load data
  new.dat.small <- make_newdat(covar.data = covar.data, sizeLevel = "small")
  new.dat.large <- make_newdat(covar.data = covar.data, sizeLevel = "large")
  
  # plotting params
  
  # ... for RCPs
  col <- rcp_colors()
  lty <- rep(1, 5)
  group_names <- paste0("Group ",1:5)
  
  # ... for variables
  variable <- c("Ca","Fe","Mn","Zn")
  xlab <- list("Calcium (%)","Iron (%)","Manganese (%)","Zinc (%)")
  
  # ... for sizeLevel
  new.dat.list <- list(small = new.dat.small, large = new.dat.large)
  
  # ... for panels
  ylab <- "Membership probability"
  npanels <- length(new.dat.list) * length(variable)
  title.mat <- matrix(letters[(1+6):(npanels+6)], nrow = length(new.dat.list), ncol = length(variable), byrow=T)
  legend.mat <- xaxis.mat <- yaxis.mat <- ymtext1.mat <- ymtext2.mat <- matrix(FALSE, nrow = length(new.dat.list), ncol = length(variable))
  xaxis.mat[2,] <- TRUE
  yaxis.mat[,1] <- TRUE
  ymtext1.mat[1,1] <- TRUE
  ymtext2.mat[2,1] <- TRUE
  
  # save
  pdf(file = "output/rcp_results/partialEffects_theRest2.pdf", width = 7.5, height = 5)
  par(mfcol=c(2,4), cex=.8) # tcl = -0.3
  par(oma=c(5,6,1,1)) # outer margins: bottom, left, top, right
  par(mar=c(0,0,1,0)+.5) # panel margins bottom, left, top, right
  
  for (i in 1:length(variable)) {
    
    for(k in 1:length(new.dat.list)){
      
      title <- title.mat[k,i]
      legend <- legend.mat[k,i]
      xaxis <- xaxis.mat[k,i]
      yaxis <- yaxis.mat[k,i]
      ymtext1 <- ymtext1.mat[k,i]
      ymtext2 <- ymtext2.mat[k,i]
      
      ParPlotMany.regimix(variable=variable[i], xlab=xlab[[i]],
                          fm=fm_rcp, 
                          newdata=new.dat.list[[k]],
                          covar.data=covar.data, 
                          RCPs=1:5, group_names=group_names, col=col, lty=lty, 
                          ylab=ylab, title=title, 
                          legend=legend, xaxis=xaxis, yaxis=yaxis, 
                          ymtext1=ymtext1, ymtext2=ymtext2)
      
    }
  }
  dev.off()
  
}


#---------------------------------------------------------#
# characterize RCP taxa

plot_mostabundTaxa_perRCP <- function(fm_rcp, taxAndFunguild){
  
  # the expected species abundances are stored in the fm$mus component of the fitted model
  # a little tricky to get - a hang over from the fit process is that the array is much much bigger than the information it contains
  # all the sites have contant species profile for each RCP, so you only need the first row of each
  
  # transpose so its RCP x species
  mus <- data.frame(t(fm_rcp$mus[1,,]))
  names(mus) <- fm_rcp$names$spp
  
  # reformat so it's easier to work with and add taxon info
  mus.df <- data.frame(t(mus))
  colnames(mus.df) <- paste0("RCP", 1:5)
  mus.df <- data.frame(OTUId = row.names(mus.df), mus.df)
  tax.indx <- unique(taxAndFunguild[,c("OTUId","OTUId_ann","phylum","family","Trophic.Mode", "Guild", "genus")])
  mus.df %>%
    gather(key = "RCP", value = "expectedAbund", -OTUId) %>%
    left_join(tax.indx) -> mus.df.l
  
  # rename RCP groups
  mus.df.l$RCP<- recode(mus.df.l$RCP, `RCP1`="Group 1", `RCP2`="Group 2", `RCP3`="Group 3", `RCP4`="Group 4", `RCP5`="Group 5")
  
  # remove values where expected abundance is below 150 reads
  mus.df.l %>%
    filter(expectedAbund > 100) %>%
    mutate(label = ifelse(grepl("OTU", OTUId_ann) == TRUE, "", OTUId_ann)) -> plot.df
  
  # trophic mode colors
  color.vec <- trophic.mode_colors()
  mytheme <- make_ggplot_theme()
  p <- ggplot(plot.df, aes(x = RCP, y = log10(expectedAbund), 
                           shape = phylum, color = Trophic.Mode)) +
    geom_point(alpha = 0.8) + 
    geom_text(aes(label = label), hjust = 0, nudge_x = 0.05, size = 2.5) + 
    expand_limits(x = 6) +
    scale_x_discrete(expand = c(0,.2)) +
    scale_shape_discrete(name = "Phylum") +
    scale_color_manual(name = "Trophic mode", values = color.vec) +
    xlab("OTU group") + ylab("Expected OTU abundance (log-transformed)") +
    mytheme + theme(legend.position="top", legend.box = "vertical", 
                    legend.key.height =unit(0,"line")) +
    guides(color = guide_legend(order=2),
           shape = guide_legend(order=1))
  p
  
  
  ggsave(file = "output/rcp_results/expectedAbund_perRCP.pdf", width = 8, height = 8)
  
  # #by Trophic.Mode
  # ggplot(subset(plot.df, Trophic.Mode %in% c("Pathotroph","Saprotroph")), 
  #               aes(x = RCP, y = log(expectedAbund), color = Trophic.Mode)) +
  #   geom_jitter(width=.2)
  
}

table_mostabundTaxa_perRCP <- function(fm_rcp, taxAndFunguild){
  
  # the expected species abundances are stored in the fm$mus component of the fitted model
  # a little tricky to get - a hang over from the fit process is that the array is much much bigger than the information it contains
  # all the sites have contant species profile for each RCP, so you only need the first row of each
  
  # transpose so its RCP x species
  mus <- data.frame(t(fm_rcp$mus[1,,]))
  names(mus) <- fm_rcp$names$spp
  
  # # and the expected cluster membership for each site
  # post_probs <- data.frame(round(fm_rcp$postProbs, 3))
  # post_probs
  # names(post_probs) <- fm_rcp$names$RCPs
  
  # reformat so it's easier to work with and add taxon info
  mus.df <- data.frame(t(mus))
  colnames(mus.df) <- paste0("RCP", 1:5)
  mus.df <- data.frame(OTUId = row.names(mus.df), mus.df)
  tax.indx <- unique(taxAndFunguild[,c("OTUId","OTUId_ann","phylum","Trophic.Mode", "Guild")])
  mus.df %>%
    gather(key = "RCP", value = "expectedAbund", -OTUId) %>%
    left_join(tax.indx) -> mus.df.l
  
  # make summary table
  mus.df.l %>%
    select(RCP, OTUId_ann, phylum, Trophic.Mode, expectedAbund) %>%
    filter(RCP == "RCP1") %>%
    arrange(desc(expectedAbund)) %>%
    slice(1:5) -> top5.rcp1
  
  mus.df.l %>%
    select(RCP, OTUId_ann, phylum, Trophic.Mode, expectedAbund) %>%
    filter(RCP == "RCP2") %>%
    arrange(desc(expectedAbund)) %>%
    slice(1:5) -> top5.rcp2
  
  mus.df.l %>%
    select(RCP, OTUId_ann, phylum, Trophic.Mode, expectedAbund) %>%
    filter(RCP == "RCP3") %>%
    arrange(desc(expectedAbund)) %>%
    slice(1:5) -> top5.rcp3
  
  mus.df.l %>%
    select(RCP, OTUId_ann, phylum, Trophic.Mode, expectedAbund) %>%
    filter(RCP == "RCP4") %>%
    arrange(desc(expectedAbund)) %>%
    slice(1:5) -> top5.rcp4
  
  mus.df.l %>%
    select(RCP, OTUId_ann, phylum, Trophic.Mode, expectedAbund) %>%
    filter(RCP == "RCP5") %>%
    arrange(desc(expectedAbund)) %>%
    slice(1:5) -> top5.rcp5
  
  top5.df <- rbind(top5.rcp1, top5.rcp2, top5.rcp3, top5.rcp4, top5.rcp5)
  
  write.csv(top5.df, file = "output/rcp_results/top5_perRCP.csv", row.names = FALSE)
  
}

