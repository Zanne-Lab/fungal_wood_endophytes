# boral - role of wood traits on fungal OTU abundances
#######


summary_Xcoefs <- function(Xcoefs.df){
  
  #count the number of signif positive and negative coefficents by run and term
  Xcoefs.df %>%
    mutate(pos = ifelse(signif == "yes" & X > 0,
                        "yes", "no")) %>%
    mutate(neg = ifelse(signif == "yes" & X < 0,
                        "yes", "no")) %>%
    group_by(term) %>%
    summarize(num_total = length(signif),
              num_signif = sum(signif == "yes"),
              num_pos = sum(pos == "yes"),
              num_neg = sum(neg == "yes")) %>%
    mutate(percPos = round((num_pos / num_total) * 100, digits = 2)) %>%
    mutate(percNeg = round((num_neg / num_total) * 100, digits = 2)) %>%
    mutate(percNonSignif = round(((num_total - num_signif) / num_total) * 100, digits = 2)) -> summary.Xcoefs.df

  return(summary.Xcoefs.df)
  
}

notResponsive_OTUs <- function(Xcoefs.df){
  
  Xcoefs.df %>%
    mutate(signif.tf = ifelse(signif == "yes", TRUE, FALSE)) %>%
    select(OTUId, term, signif.tf) %>%
    spread(key = term, value = signif.tf) -> df
  rownames(df)<-df$OTUId
  numSignif<-rowSums(df[,-1])
  num.nr <- length(numSignif[numSignif == 0])
  
  return(num.nr)
  }

write_summary_Xcoefs <- function(fit.list, allXs){
  
  if(allXs == TRUE){
    modelID2 <- "allX"
  }else{
    modelID2 <- "selectX"
  }
  
  Xcoefs.list <- lapply(fit.list, function(x) x$Xcoefs.df)
  Xcoefs.summary.list <- lapply(Xcoefs.list, summary_Xcoefs)
  Xcoefs.summary.df <- list_to_df(Xcoefs.summary.list)
  
  fileName <- paste0("output/boral_roleOfTraits/summary_Xcoefs_",modelID2, ".csv")
  write.csv(Xcoefs.summary.df, file=fileName, row.names = FALSE)
  
  
  # calculate what % of OTUs do not respond significantly to even 1 X var
  num.nr.list <- lapply(Xcoefs.list, notResponsive_OTUs)
  num.nr.list
}

summary_OTUID_Xcoefs <- function(Xcoefs.df, taxAndFunguild){
  
  #make OTU look up tab
  taxAndFunguild %>%
    select(OTUId, OTUId_ann, genus, phylum, Trophic.Mode) -> tax.indx
  
  #make summary table
  Xcoefs.df %>%
    left_join(tax.indx) %>%
    filter(signif == "yes") %>%
    filter(genus != "unclassified") %>%
    mutate(X.round = round(X, digits = 4)) %>%
    select(OTUId_ann, term, X.round, phylum, genus, Trophic.Mode) %>%
    spread(key = term, value = X.round) %>%
    rename("Phylum" = "phylum",
           "Genus" = "genus",
           "Guild" = "Trophic.Mode",
           "BT" = "barkthick",
           "WP" = "waterperc",
           "Size" = "sizesmall") -> Xcoef.tab
  
  # density column dropped sometimes, so force it back in if missing
  test <- sum(colnames(Xcoef.tab) %in% "density")
  if(test == 0){
    Xcoef.tab$density <- NA
    Xcoef.tab %>%
      select(OTUId_ann, Phylum, Genus, Guild, BT, C, Ca, density, Fe, K, Mn, N, P, Size, WP, Zn) -> Xcoef.tab
  }

  #shorten codes and arrange rows
  # Xcoef.tab$Phylum <- recode(Xcoef.tab$Phylum, 
  #                            "Ascomycota"="A", 
  #                            "Basidiomycota"="B")
  # Xcoef.tab$Guild <- recode(Xcoef.tab$Guild, 
  #                           "Saprotroph" = "Sapr",
  #                           "Pathotroph" = "Path",
  #                           "Symbiotroph" = "Symb",
  #                           "Pathotroph-Symbiotroph" = "Path-Symb",
  #                           "Pathotroph-Saprotroph" = "Path-Sapr")
  # Xcoef.tab %>% mutate(Guild = ifelse(Guild == "unclassified",
  #                                     NA, Guild)) %>%
  #   group_by(Phylum, Guild, Genus) %>% arrange(Size, .by_group = TRUE) -> Xcoef.tab
  
  return(Xcoef.tab)

}

summary_OTUID_nonsignif_Xcoefs <- function(Xcoefs.df, taxAndFunguild){
  
  #make OTU look up tab
  taxAndFunguild %>%
    select(OTUId, OTUId_ann, genus, phylum, Trophic.Mode) -> tax.indx
  
  #make summary table
  Xcoefs.df %>%
    left_join(tax.indx) %>%
    filter(genus != "unclassified") %>%
    mutate(X.round = round(X, digits = 4)) %>%
    select(OTUId, term, X.round) %>%
    spread(key = term, value = X.round) %>%
    rename("BT" = "barkthick",
           "WP" = "waterperc",
           "Size" = "sizesmall") -> Xcoef.tab
  
  #select(OTUId_ann, phylum, genus, Trophic.Mode, sizesmall, density, barkthick, waterperc, C, N, P, Ca, Fe, K, Mn, Zn) %>%
  
  #shorten codes and arrange rows
  # Xcoef.tab$Phylum <- recode(Xcoef.tab$Phylum, 
  #                            "Ascomycota"="A", 
  #                            "Basidiomycota"="B")
  # Xcoef.tab$Guild <- recode(Xcoef.tab$Guild, 
  #                           "Saprotroph" = "Sapr",
  #                           "Pathotroph" = "Path",
  #                           "Symbiotroph" = "Symb",
  #                           "Pathotroph-Symbiotroph" = "Path-Symb",
  #                           "Pathotroph-Saprotroph" = "Path-Sapr")
  # Xcoef.tab %>% mutate(Guild = ifelse(Guild == "unclassified",
  #                                     NA, Guild)) %>%
  #   group_by(Phylum, Guild, Genus) %>% arrange(Size, .by_group = TRUE) -> Xcoef.tab
  
  return(Xcoef.tab)
  
}

# signif result must be in 9/12 runs or more
plot_summary_Xcoefs_byOTUId <- function(fit.list, taxAndFunguild, allXs){
  
  if(allXs == TRUE){
    modelID2 <- "allX"
    xvar.order <- c("Size","density","WP","BT","C","N","P","Ca","Fe","K","Mn","Zn")
    fig.width <- 10
  }else{
    modelID2 <- "selectX"
    xvar.order <- c("Size","WP","BT","C")
    fig.width <- 5
  }
  
  Xcoefs.list <- lapply(fit.list, function(x) x$Xcoefs.df)
  Xcoefs.summary.list <- lapply(Xcoefs.list, function(x){
    result<- summary_OTUID_Xcoefs(Xcoefs.df = x, taxAndFunguild = taxAndFunguild)
    result.df<-data.frame(result)
    return(result.df)
    })
  Xcoefs.summary.df <- list_to_df(Xcoefs.summary.list)
  
  # create summary plot
  Xcoefs.summary.df %>%
    gather(key = "term", value = "estimate", -c(OTUId_ann, Phylum, Genus, Guild, source), na.rm = T) -> Xcoefs.summary.df.l
  
  # identify OTU and Xvar estimates that flip sign and remove from the dataset
  Xcoefs.summary.df.l %>%
    mutate(estimate_sign = ifelse(estimate > 0, "positive","negative")) %>%
    group_by(OTUId_ann, term) %>%
    summarize(numSign = length(unique(estimate_sign))) %>%
    filter(numSign != 1) -> flip.signif.indx
  
  Xcoefs.summary.df.l %>%
    left_join(flip.signif.indx) %>%
    filter(is.na(numSign)) -> Xcoefs.summary.df.l.noflips
  
  # identify signif OTU and Xvar estimate that only show up in some of the runs
  Xcoefs.summary.df.l.noflips %>%
    group_by(OTUId_ann, term) %>%
    summarize(numRuns = length(unique(source))) -> all.runs.indx
  Xcoefs.summary.df.l.noflips %>%
    left_join(all.runs.indx) %>%
    filter(numRuns >= 9) -> Xcoefs.summary.df.l.noflips.allruns
  
  # order the terms
  df <- Xcoefs.summary.df.l.noflips.allruns
  df$term <- factor(df$term, levels = xvar.order)
  
  # summarize across runs
  df %>%
    group_by(OTUId_ann, term) %>%
    summarize(meanEst = mean(estimate),
              seEst = sd(estimate)/sqrt(length(estimate))) -> summ.df
  indx <- unique(Xcoefs.summary.df[,c("OTUId_ann","Phylum","Genus","Guild")])
  summ.df %>%
    left_join(indx) -> summ.df.ann
  
  # order the OTUs
  if(allXs == TRUE){
    summ.df$term <- recode(summ.df$term, "density"="Density", "WP"="Water", "BT"="Bark")
    summ.df %>%
      select(OTUId_ann, term, meanEst) %>%
      spread(key = term, value = meanEst) %>%
      arrange(Size, Density, Water, Bark, C, N, P, Ca, Fe, K, Mn, Zn) -> otu.order
    summ.df.ann$OTUId_ann <- factor(summ.df.ann$OTUId_ann, levels = rev(otu.order$OTUId_ann)) 
    summ.df.ann$term <- recode(summ.df.ann$term, "density"="Density", "WP"="Water", "BT"="Bark")
  }else{
    summ.df %>%
      select(OTUId_ann, term, meanEst) %>%
      spread(key = term, value = meanEst) %>%
      arrange(Size, WP) -> otu.order
    summ.df.ann$OTUId_ann <- factor(summ.df.ann$OTUId_ann, levels = rev(otu.order$OTUId_ann)) 
  }
  
  color.vec<- trophic.mode_colors()
  mytheme <- make_ggplot_theme()
  ggplot(summ.df.ann, aes(y = OTUId_ann, x = meanEst, shape = Phylum, color = Guild)) +
    geom_point(size = 3) +
    geom_errorbarh(aes(xmin = meanEst - seEst, xmax = meanEst + seEst)) +
    facet_wrap(~term, scales = "free_x", nrow = 1) +
    scale_color_manual(name = "Trophic mode", values = color.vec) +
    geom_vline(aes(xintercept = 0), linetype = 2) +
    mytheme + theme(panel.grid.major = element_line("grey70", size = 0.2)) +
    theme(legend.position="top", legend.box = "vertical",
          legend.key.height=unit(0,"line")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    xlab("Coefficient estimate") + ylab("OTU") +
    guides(color = guide_legend(order=2),
           shape = guide_legend(order=1))
  
  fileName <- paste0("output/boral_roleOfTraits/signifXcoefs_byOTUId_", modelID2, ".pdf")
  ggsave(file = fileName, width = fig.width, height = 6)

}
