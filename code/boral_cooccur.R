# boral - cooccurance results
#######

#---------------------------------------------------------#
# prep data

annotate_cor_withOTUInfo <- function(cor.df, taxAndFunguild){
  
  # taxa index
  taxAndFunguild %>%
    select(OTUId, genus, phylum, Trophic.Mode) -> tax.indx
  
  # link to otu1
  tax.indx %>%
    rename("otu1"="OTUId") -> tax.indx
  cor.df %>%
    left_join(tax.indx) %>%
    rename("genus_1"="genus",
           "phylum_1"="phylum",
           "troph_1"="Trophic.Mode") %>%
    mutate(genus_1_unclass = ifelse(genus_1 == "unclassified",
                                    "yes","no")) -> cor.df.1ann
  
  #link to otu2
  tax.indx %>%
    rename("otu2"="otu1") -> tax.indx
  cor.df.1ann %>%
    left_join(tax.indx) %>%
    rename("genus_2"="genus",
           "phylum_2"="phylum",
           "troph_2"="Trophic.Mode") %>%
    mutate(genus_2_unclass = ifelse(genus_2 == "unclassified",
                                    "yes","no")) -> cor.df.ann
  
  #classify pairs
  cor.df.ann %>%
    mutate(pair_haveGenusIDs = ifelse(genus_1_unclass == "no" & genus_2_unclass == "no",
                                      "yes", "no")) %>%
    mutate(pair_samePhylum = ifelse(phylum_1 == phylum_2,
                                    "yes","no")) %>%
    mutate(pair_sameTroph = ifelse(troph_1 == troph_2,
                                   "yes","no")) -> cor.df.annpairs
  
  #make same phylum more specific
  cor.df.annpairs %>%
    mutate(pair_samePhylum = ifelse(!phylum_1 %in% c("Ascomycota","Basidiomycota") | !phylum_2 %in% c("Ascomycota","Basidiomycota"),
                                    "notAsco-Basidio", pair_samePhylum)) -> cor.df.annpairs
  
  #make same troph more specific
  usethese.trophs <- c("Pathotroph","Saprotroph")
  cor.df.annpairs %>%
    mutate(pair_sameTroph = ifelse(!troph_1 %in% usethese.trophs | !troph_2 %in% usethese.trophs,
                                   "notPatho-Sapro", pair_sameTroph)) -> cor.df.annpairs
  
  #select cols
  cor.df.annpairs %>%
    select(source, otu1, otu2, 
           enviro_cor, enviro_signif, residual_cor, residual_signif,
           phylum_1, genus_1, phylum_2, genus_2,
           pair_haveGenusIDs, pair_samePhylum, pair_sameTroph) -> cor.df.annpairs
  
  return(cor.df.annpairs)
  
}

annotate_cor_withOTUInfo_v <- function(cor.df, taxAndFunguild){
  
  # taxa index
  taxAndFunguild %>%
    select(OTUId, OTUId_ann, genus, phylum, Trophic.Mode) -> tax.indx
  
  # link to otu1
  tax.indx %>%
    rename("otu1"="OTUId") -> tax.indx
  cor.df %>%
    left_join(tax.indx) %>%
    rename("OTUIdann_1"="OTUId_ann",
           "genus_1"="genus",
           "phylum_1"="phylum",
           "troph_1"="Trophic.Mode") %>%
    mutate(genus_1_unclass = ifelse(genus_1 == "unclassified",
                                    "yes","no")) -> cor.df.1ann
  
  #link to otu2
  tax.indx %>%
    rename("otu2"="otu1") -> tax.indx
  cor.df.1ann %>%
    left_join(tax.indx) %>%
    rename("OTUIdann_2"="OTUId_ann",
           "genus_2"="genus",
           "phylum_2"="phylum",
           "troph_2"="Trophic.Mode") %>%
    mutate(genus_2_unclass = ifelse(genus_2 == "unclassified",
                                    "yes","no")) -> cor.df.ann
  
  #classify pairs
  cor.df.ann %>%
    mutate(pair_haveGenusIDs = ifelse(genus_1_unclass == "no" & genus_2_unclass == "no",
                                      "yes", "no")) %>%
    mutate(pair_samePhylum = ifelse(phylum_1 == phylum_2,
                                    "yes","no")) %>%
    mutate(pair_sameTroph = ifelse(troph_1 == troph_2,
                                   "yes","no")) -> cor.df.annpairs
  
  #make same phylum more specific
  cor.df.annpairs %>%
    mutate(pair_samePhylum = ifelse(!phylum_1 %in% c("Ascomycota","Basidiomycota") | !phylum_2 %in% c("Ascomycota","Basidiomycota"),
                                    "notAsco-Basidio", pair_samePhylum)) -> cor.df.annpairs
  
  #make same troph more specific
  usethese.trophs <- c("Pathotroph","Saprotroph")
  cor.df.annpairs %>%
    mutate(pair_sameTroph = ifelse(!troph_1 %in% usethese.trophs | !troph_2 %in% usethese.trophs,
                                   "notPatho-Sapro", pair_sameTroph)) -> cor.df.annpairs

  #prune columns
  cor.df.annpairs %>%
    select(otu1, OTUIdann_1, otu2, OTUIdann_2, enviro_cor, residual_cor, source, phylum_1, troph_1, phylum_2, troph_2, 
           pair_haveGenusIDs, enviro_signif, residual_signif) -> result
  
  return(result)
  
}


#---------------------------------------------------------#
# How different are correlations across runs?

enviro.cor.btwRuns <- function(fit.list, taxAndFunguild, allXs){
  
  if(allXs == TRUE){
    modelID2 <- "allX"
  }else{
    modelID2 <- "selectX"
  }
  
  cor.list <- lapply(fit.list, function(x) x$cor.df)
  cor.df <- list_to_df(cor.list)
  cor.df.ann <- annotate_cor_withOTUInfo(cor.df = cor.df, taxAndFunguild = taxAndFunguild)
  cor.df.ann %>%
    mutate(otupair = paste0(otu1, otu2)) -> cor.df.ann
  
  # how many pairs in total?
  numOTUpairs <- length(unique(cor.df.ann$otupair))
  
  cor.df.ann %>%
    select(source, otupair, enviro_cor, enviro_signif) %>%
    group_by(otupair) %>%
    summarize(numSignifLevels = length(unique(enviro_signif)),
              signifLevels = paste(unique(enviro_signif), collapse = "_")) -> summ.enviro
  
  # look at OTUs with more than 1 signif category
  summ.enviro %>%
    filter(numSignifLevels !=1) -> tmp
  
  # look at OTUs that classified as signif positive and negative
  summ.enviro %>%
    filter(grepl("positive", signifLevels) & grepl("negative", signifLevels)) -> tmp
  numDiff <-dim(tmp)[1]
  percDiff <- (numDiff / numOTUpairs) *100
  
  # # calculate the difference between runs
  # cor.df.ann %>%
  #   select(source, otupair, enviro_cor) %>%
  #   filter(otupair %in% tmp$otupair) %>%
  #   group_by(otupair) %>%
  #   summarize(min.cor = min(enviro_cor),
  #             max.cor = max(enviro_cor),
  #             diff = abs(max.cor - min.cor)) %>%
  #   select(otupair, diff) -> diff.indx
  # 
  # # subset and annotate the data
  # cor.df.ann %>%
  #   select(source, otupair, enviro_cor, enviro_signif) %>%
  #   filter(otupair %in% tmp$otupair) %>%
  #   left_join(diff.indx) -> plot.df
  # 
  # # plot the distribution of differences
  # PAIR<-unique(plot.df$otupair)
  # title <- paste(length(PAIR), " of ", numOTUpairs, 
  #                " (", round((length(PAIR)/numOTUpairs) * 100, digits = 0), "%) ",
  #                "OTU pairs flip signif", sep ="")
  # hist.envir <- ggplot(diff.indx, aes(x=diff)) +
  #   geom_histogram() + xlab("Abs. diff. in correlation values between runs") +
  #   ggtitle(title) +
  #   geom_vline(aes(xintercept = 0.5), color = 2)
  # #hist.envir
  # 
  # # plot the correlation values for OTU pairs where the difference between runs is greater that 0.5
  # plot.df -> plot.df.sub
  # badones.envir <- ggplot(plot.df.sub, aes(y = reorder(otupair, diff), x = enviro_cor)) +
  #   geom_point(aes(color = enviro_signif)) +
  #   geom_line() +
  #   geom_vline(aes(xintercept = 0), linetype = 2) +
  #   guides(color = F) +
  #   xlab("Correlation") + ylab("OTU pair")
  # #badones.envir
  
  #fileName <- paste0("output/boral_cooccur/runVariation_enviro_", modelID2, ".pdf")
  #pdf(fileName, width = 12, height = 8)
  #grid.arrange(hist.envir, badones.envir, ncol=2)
  #dev.off()
  
  return(percDiff) # percent of OTUpairs that flip significant signs (e.g. some runs they are positive, others they are negative)
  
}

residual.cor.btwRuns <- function(fit.list, taxAndFunguild, allXs){
  
  if(allXs == TRUE){
    modelID2 <- "allX"
  }else{
    modelID2 <- "selectX"
  }
  
  cor.list <- lapply(fit.list, function(x) x$cor.df)
  cor.df <- list_to_df(cor.list)
  cor.df.ann <- annotate_cor_withOTUInfo(cor.df = cor.df, taxAndFunguild = taxAndFunguild)
  cor.df.ann %>%
    mutate(otupair = paste0(otu1, otu2)) -> cor.df.ann
  
  # how many pairs in total?
  numOTUpairs <- length(unique(cor.df.ann$otupair))
  
  cor.df.ann %>%
    select(source, otupair, residual_cor, residual_signif) %>%
    group_by(otupair) %>%
    summarize(numSignifLevels = length(unique(residual_signif)),
              signifLevels = paste(unique(residual_signif), collapse = "_")) -> summ.residual
  
  # look at OTUs with more than 1 signif category
  summ.residual %>%
    filter(numSignifLevels !=1) -> tmp
  
  # look at OTUs that classified as signif positive and negative
  summ.residual %>%
    filter(grepl("positive", signifLevels) & grepl("negative", signifLevels)) -> tmp
  numDiff <-dim(tmp)[1]
  percDiff <- (numDiff / numOTUpairs) *100
  percDiff
  
  # # calculate the difference between runs
  # cor.df.ann %>%
  #   select(source, otupair, residual_cor) %>%
  #   filter(otupair %in% tmp$otupair) %>%
  #   group_by(otupair) %>%
  #   summarize(min.cor = min(residual_cor),
  #             max.cor = max(residual_cor),
  #             diff = abs(max.cor - min.cor)) %>%
  #   select(otupair, diff) -> diff.indx
  # 
  # # subset and annotate the data
  # cor.df.ann %>%
  #   select(source, otupair, residual_cor, residual_signif) %>%
  #   filter(otupair %in% tmp$otupair) %>%
  #   left_join(diff.indx) -> plot.df
  # 
  # # plot the distribution of differences
  # PAIR<-unique(plot.df$otupair)
  # title <- paste(length(PAIR), " of ", numOTUpairs, 
  #                " (", round((length(PAIR)/numOTUpairs) * 100, digits = 0), "%) ",
  #                "OTU pairs flip signif", sep ="")
  # hist.residual <- ggplot(diff.indx, aes(x=diff)) +
  #   geom_histogram() + xlab("Abs. diff. in correlation values between runs") +
  #   ggtitle(title) +
  #   geom_vline(aes(xintercept = 0.15), color = 2)
  # #hist.residual
  # 
  # # plot the correlation values for OTU pairs where the difference between runs is greater that 0.5
  # plot.df -> plot.df.sub
  # badones.residual <- ggplot(plot.df.sub, aes(y = reorder(otupair, diff), x = residual_cor)) +
  #   geom_point(aes(color = residual_signif)) +
  #   geom_line() +
  #   geom_vline(aes(xintercept = 0), linetype = 2) +
  #   guides(color = F) +
  #   xlab("Correlation") + ylab("OTU pair")
  # #badones.residual
  # 
  # fileName <- paste0("output/boral_cooccur/runVariation_resdiual_", modelID2, ".pdf")
  # pdf(fileName, width = 12, height = 8)
  # grid.arrange(hist.residual, badones.residual, ncol=2)
  # dev.off()
  
  return(percDiff) # percent of OTUpairs that flip significant signs (e.g. some runs they are positive, others they are negative)
  
}


#---------------------------------------------------------#
# plot distribution of correlation values
plot_cor_distributions_allruns <- function(fit.list, taxAndFunguild, allXs){
  
  if(allXs == TRUE){
    modelID2 <- "allX"
  }else{
    modelID2 <- "selectX"
  }
  
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
  
  fileName <- paste0("output/boral_cooccur/cor_distributions_allruns_",modelID2,".pdf")
  pdf(file=fileName, width=12, height=6)
  grid.arrange(p.env, 
               p.res, 
               nrow = 2)
  dev.off()
  
}

plot_cor_distributions <- function(fit.list, taxAndFunguild, allXs){
  
  if(allXs == TRUE){
    modelID2 <- "allX"
  }else{
    modelID2 <- "selectX"
  }
  
  cor.list <- lapply(fit.list, function(x) x$cor.df)
  cor.df <- list_to_df(cor.list)
  cor.df.ann <- annotate_cor_withOTUInfo(cor.df = cor.df, taxAndFunguild = taxAndFunguild)
  
  #just use the first run's data
  cor.df.ann %>%
    filter(source == "run1") -> cor.df.ann
  
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
  cor.df.ann %>%
    group_by(enviro_signif) %>%
    summarize(num = length(enviro_signif)) -> tmp
  total <- sum(tmp[,2])
  percPos <- round((tmp[tmp$enviro_signif=="positive", "num"]/total)*100, digits = 0)
  percNeg <- round((tmp[tmp$enviro_signif=="negative", "num"]/total)*100, digits = 0)
  p.env<-ggplot(cor.df.ann, aes(x = enviro_cor, fill = enviro_signif, alpha = pair_haveGenusIDs)) + 
    geom_histogram(bins = 30) +
    geom_vline(xintercept = 0, linetype=2) +
    xlab("Shared environment correlation") + ylab("Frequency") +
    scale_fill_manual(name = "Signif.", values = colorVec) +
    scale_alpha_manual(name = "Genera IDs", values = c(0.5, 1)) +
    xlim(c(-1,1)) +
    mytheme + guides(fill=FALSE, alpha=FALSE) +
    annotate("text", y=300, x=-0.45, label=paste(percNeg, "%", sep="")) +
    annotate("text", y=300, x=0.5, label=paste(percPos ,"%", sep=""))
  p.env
  
  #residual
  cor.df.ann %>%
    group_by(residual_signif) %>%
    summarize(num = length(residual_signif)) -> tmp
  total <- sum(tmp[,2])
  percPos <- round((tmp[tmp$residual_signif=="positive", "num"]/total)*100, digits = 0)
  percNeg <- round((tmp[tmp$residual_signif=="negative", "num"]/total)*100, digits = 0)
  p.res<-ggplot(cor.df.ann, aes(x = residual_cor, fill = residual_signif, alpha = pair_haveGenusIDs)) + 
    geom_histogram(bins = 30) +
    geom_vline(xintercept = 0, linetype=2) +
    xlab("Residual correlation") + ylab("Frequency") +
    scale_fill_manual(name = "Signif.", values = colorVec) +
    scale_alpha_manual(name = "Genera IDs", values = c(0.5, 1)) +
    xlim(c(-1,1)) +
    mytheme + guides(fill=FALSE, alpha=FALSE) +
    annotate("text", y=300, x=-0.3, label=paste(percNeg,"%", sep="")) +
    annotate("text", y=300, x=0.5, label=paste(percPos,"%", sep="")) 
  p.res
    
  fileName <- paste0("output/boral_cooccur/cor_distributions_",modelID2,".pdf")
  pdf(file=fileName, width=6, height=3)
  grid.arrange(p.env, 
               p.res, 
               ncol = 2)
  dev.off()
  
}

write_summary_cor_distribution <- function(fit.list, taxAndFunguild, allXs){
  
  if(allXs == TRUE){
    modelID2 <- "allX"
  }else{
    modelID2 <- "selectX"
  }
  
  cor.list <- lapply(fit.list, function(x) x$cor.df)
  cor.df <- list_to_df(cor.list)
  cor.df.ann <- annotate_cor_withOTUInfo(cor.df = cor.df, taxAndFunguild = taxAndFunguild)
  
  #just use the first run's data
  cor.df.ann %>%
    filter(source == "run1") %>%
    select(source, otu1, otu2, enviro_cor, residual_cor) %>%
    gather(key = "key",value = "value", -c(source, otu1, otu2)) -> df

  #summarize the enviro and residual distributions
  df %>%
    group_by(key) %>%
    summarize(medval = median(value),
              meanval = mean(value),
              seval = sd(value)/sqrt(sum(!is.na(value))),
              lowerlim = range(value)[1],
              upperlim = range(value)[2]) -> df.summ
  
  fileName <- paste0("output/boral_cooccur/cor_distribution_summarystats_",modelID2,".csv")
  write.csv(df.summ, file=fileName)
  
}

#---------------------------------------------------------#
# plot frequence of significant positive and negative correlation values by group

plot_corFreq_phylo <- function(fit.list, taxAndFunguild, allXs){
  
  if(allXs == TRUE){
    modelID2 <- "allX"
  }else{
    modelID2 <- "selectX"
  }
  
  cor.list <- lapply(fit.list, function(x) x$cor.df)
  cor.df <- list_to_df(cor.list)
  cor.df.ann <- annotate_cor_withOTUInfo(cor.df = cor.df, taxAndFunguild = taxAndFunguild)
  
  #just use the first run's data
  cor.df.ann %>%
    filter(source == "run2") -> cor.df.ann
  
  #summarize frequency of pos, neg, zero correlations by phylo pair
  cor.df.ann %>%
    filter(pair_samePhylum != "notAsco-Basidio") %>%
    group_by(pair_samePhylum) %>%
    summarize(num_total = length(enviro_signif),
              num_env_pos = sum(enviro_signif == "negative"),
              num_env_neg = sum(enviro_signif == "positive"),
              num_res_pos = sum(residual_signif == "negative"),
              num_res_neg = sum(residual_signif == "positive")) %>%
    mutate(env_percPos = round((num_env_pos / num_total) * 100, digits = 2)) %>%
    mutate(env_percNeg = round((num_env_neg / num_total) * 100, digits = 2)) %>%
    mutate(res_percPos = round((num_res_pos / num_total) * 100, digits = 2)) %>%
    mutate(res_percNeg = round((num_res_neg / num_total) * 100, digits = 2)) %>%
    select(pair_samePhylum, num_total, env_percPos, env_percNeg, res_percPos, res_percNeg) %>%
    gather(key = "key", value = "percent", 3:6) %>%
    separate(key, into = c("corType","sign")) -> summ
  
  #clean sign factor
  summ$sign <- recode(summ$sign, 
                      "percNeg" = "Negative",
                      "percPos" = "Positive")
  
  #clean pair factor
  tot <- unique(summ[,c("pair_samePhylum","num_total")])
  tot
  summ$pair_samePhylum <- recode(summ$pair_samePhylum, 
                                 "no" = "Between (n=2813)",
                                 "yes" = "Within (n=5062)")
  
  #plotting params
  mytheme <- make_ggplot_theme()
  colorVec <- cor_colors()
  names(colorVec) <- c("Negative","Positive")
  
  #enviro
  summ %>%
    filter(corType == "env") -> plotdf
  p.env <- ggplot(plotdf, aes(x = sign, y = percent, fill = sign, alpha = pair_samePhylum)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.8, col = 1) +
    ylab("Frequency (%)") + xlab("Shared environment\ncorrelation sign") +
    scale_fill_manual(values = colorVec) +
    scale_alpha_manual(name = "Group", values = c(0.5, 1)) +
    mytheme + guides(fill = FALSE)
  
  #residual
  summ %>%
    filter(corType == "res") -> plotdf
  p.res <- ggplot(plotdf, aes(x = sign, y = percent, fill = sign, alpha = pair_samePhylum)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.8, col = 1) +
    ylab("") + xlab("Residual\ncorrelation sign") +
    scale_fill_manual(values = colorVec) +
    scale_alpha_manual(name = "Group", values = c(0.5, 1)) +
    mytheme + guides(fill = FALSE)
  
  #need to fix the y axis scale on these later
  fileName <- paste0("output/boral_cooccur/corFreq_phylo_",modelID2,".pdf")
  pdf(file=fileName, width=5, height=2.5)
  grid.arrange(p.env + ggtitle("c") + guides(alpha = FALSE), 
               p.res + ggtitle("d") + guides(alpha = FALSE),
               g_legend(p.res),
               ncol = 3)
  dev.off()
  
  }

plot_corFreq_troph <- function(fit.list, taxAndFunguild, allXs){
  
  if(allXs == TRUE){
    modelID2 <- "allX"
  }else{
    modelID2 <- "selectX"
  }
  
  cor.list <- lapply(fit.list, function(x) x$cor.df)
  cor.df <- list_to_df(cor.list)
  cor.df.ann <- annotate_cor_withOTUInfo(cor.df = cor.df, taxAndFunguild = taxAndFunguild)
  
  #just use the first run's data
  cor.df.ann %>%
    filter(source == "run2") -> cor.df.ann
  
  #summarize frequency of pos, neg, zero correlations by phylo pair
  cor.df.ann %>%
    filter(pair_sameTroph != "notPatho-Sapro") %>%
    group_by(pair_sameTroph) %>%
    summarize(num_total = length(enviro_signif),
              num_env_pos = sum(enviro_signif == "negative"),
              num_env_neg = sum(enviro_signif == "positive"),
              num_res_pos = sum(residual_signif == "negative"),
              num_res_neg = sum(residual_signif == "positive")) %>%
    mutate(env_percPos = round((num_env_pos / num_total) * 100, digits = 2)) %>%
    mutate(env_percNeg = round((num_env_neg / num_total) * 100, digits = 2)) %>%
    mutate(res_percPos = round((num_res_pos / num_total) * 100, digits = 2)) %>%
    mutate(res_percNeg = round((num_res_neg / num_total) * 100, digits = 2)) %>%
    select(pair_sameTroph, num_total, env_percPos, env_percNeg, res_percPos, res_percNeg) %>%
    gather(key = "key", value = "percent", 3:6) %>%
    separate(key, into = c("corType","sign")) -> summ
  
  #clean sign factor
  summ$sign <- recode(summ$sign, 
                      "percNeg" = "Negative",
                      "percPos" = "Positive")
  
  #clean pair factor
  tot <- unique(summ[,c("pair_sameTroph","num_total")])
  tot
  summ$pair_sameTroph <- recode(summ$pair_sameTroph, 
                                 "no" = "Between (n=221)",
                                 "yes" = "Within (n=214)")
  
  #plotting params
  mytheme <- make_ggplot_theme()
  colorVec <- cor_colors()
  names(colorVec) <- c("Negative","Positive")
  
  #enviro
  summ %>%
    filter(corType == "env") -> plotdf
  p.env <- ggplot(plotdf, aes(x = sign, y = percent, fill = sign, alpha = pair_sameTroph)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.8, col = 1) +
    ylab("Frequency (%)") + xlab("Shared environment\ncorrelation sign") +
    scale_fill_manual(values = colorVec) +
    scale_alpha_manual(name = "Group", values = c(0.5, 1)) +
    mytheme + guides(fill = FALSE)
  
  #residual
  summ %>%
    filter(corType == "res") -> plotdf
  p.res <- ggplot(plotdf, aes(x = sign, y = percent, fill = sign, alpha = pair_sameTroph)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.8, col = 1) +
    ylab("") + xlab("Residual\ncorrelation sign") +
    scale_fill_manual(values = colorVec) +
    scale_alpha_manual(name = "Group", values = c(0.5, 1)) +
    mytheme + guides(fill = FALSE)
  
  #need to fix the y axis scale on these later
  fileName <- paste0("output/boral_cooccur/corFreq_troph_",modelID2,".pdf")
  pdf(file=fileName, width=5, height=2.5)
  grid.arrange(p.env + ggtitle("a") + guides(alpha = FALSE), 
               p.res + ggtitle("b") + guides(alpha = FALSE),
               g_legend(p.res),
               ncol = 3)
  dev.off()

}


#---------------------------------------------------------#
# chord diagrams

make_chordDiagrams <- function(fit.list, taxAndFunguild, allXs){
  
  if(allXs == TRUE){
    modelID2 <- "allX"
  }else{
    modelID2 <- "selectX"
  }
  
  require(circlize)
  
  cor.list <- lapply(fit.list, function(x) x$cor.df)
  cor.df <- list_to_df(cor.list)
  cor.df.ann <- annotate_cor_withOTUInfo(cor.df = cor.df, taxAndFunguild = taxAndFunguild)
  cor.df.ann %>%
    mutate(otupair = paste0(otu1, otu2)) -> cor.df.ann
  
  #only use OTU pairs where both are IDd to genus
  cor.df.ann %>%
    filter(pair_haveGenusIDs == "yes") -> cor.df.ann
  
  #plotting params
  colorVec <- cor_colors()
  names(colorVec) <- c("negative","positive","none")
  
  ########################
  #enviro
  
  #identify correlations between OTU pairs that are stable
  cor.df.ann %>%
    select(source, otupair, enviro_cor, enviro_signif) %>%
    group_by(otupair) %>%
    summarize(numSignifLevels = length(unique(enviro_signif)),
              signifLevels = paste(unique(enviro_signif), collapse = "_")) -> summ.enviro
  summ.enviro$signifLevels <- recode(summ.enviro$signifLevels, negative_no = "no_negative", positive_no = "no_positive",
                                     negative_no_positive = "no_negative_positive", positive_no_negative = "no_negative_positive")
  summ.enviro %>%
    filter(signifLevels %in% c("negative","positive")) -> stable.otupairs
  # make the plotting dataframe
  cor.df.ann %>%
    filter(source == "run1") %>%
    filter(otupair %in% stable.otupairs$otupair) %>%
    select(genus_1, genus_2, enviro_cor) %>%
    mutate(col = ifelse(enviro_cor > 0, 
                        colorVec['positive'], colorVec['negative'])) %>%
    arrange(desc(enviro_cor)) -> env.frame
  dim(env.frame)
  
  ########################
  #residual
  
  #identify correlations between OTU pairs that are stable
  cor.df.ann %>%
    select(source, otupair, residual_cor, residual_signif) %>%
    group_by(otupair) %>%
    summarize(numSignifLevels = length(unique(residual_signif)),
              signifLevels = paste(unique(residual_signif), collapse = "_")) -> summ.residual
  summ.residual$signifLevels <- recode(summ.residual$signifLevels, negative_no = "no_negative", positive_no = "no_positive",
                                     negative_no_positive = "no_negative_positive", positive_no_negative = "no_negative_positive")
  summ.residual %>%
    filter(signifLevels %in% c("negative","positive")) -> stable.otupairs
  # make the plotting dataframe
  cor.df.ann %>%
    filter(source == "run1") %>%
    filter(otupair %in% stable.otupairs$otupair) %>%
    select(genus_1, genus_2, residual_cor) %>%
    mutate(col = ifelse(residual_cor > 0, 
                        colorVec['positive'], colorVec['negative'])) %>%
    arrange(desc(residual_cor)) -> res.frame
  dim(res.frame)
  
  
  ########################
  # plot
  
  #make plots
  fileName <- paste0("output/boral_cooccur/chordDiagrams_", modelID2, ".pdf")
  pdf(file=fileName, width=10, height=5, colormodel = "cmyk")
  par(mfrow=c(1,2))
  circos.par(track.height=.5)
  
  ### env
  chordDiagram(env.frame[,1:3], transparency = .2, 
               col = env.frame$col, 
               annotationTrack = c("grid"), 
               annotationTrackHeight = c(0.01,0.01),
               preAllocateTracks = 1,
               grid.col= 1, reduce = 0)
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.text(mean(xlim), ylim[1] + .001, sector.name, 
                facing = "clockwise", 
                niceFacing = TRUE, adj = c(0, 0.5), cex=.7)
  }, bg.border = NA)
  
  ### res
  chordDiagram(res.frame[,1:3], transparency = .2,
               col = res.frame$col, 
               annotationTrack = c("grid"), 
               annotationTrackHeight = c(0.01,0.01),
               preAllocateTracks = 1,
               grid.col= 1, reduce = 0)
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.text(mean(xlim), ylim[1] + .001, sector.name, 
                facing = "clockwise", 
                niceFacing = TRUE, adj = c(0, 0.5), cex=.7)
  }, bg.border = NA)
  
  dev.off()
  
  
}

investigate_enviro_chordDiagram <- function(fit.list, fit.list2, taxAndFunguild, allXs = TRUE){
  
  cor.list <- lapply(fit.list, function(x) x$cor.df)
  cor.df <- list_to_df(cor.list)
  cor.df.ann <- annotate_cor_withOTUInfo(cor.df = cor.df, taxAndFunguild = taxAndFunguild)
  cor.df.ann %>%
    mutate(otupair = paste0(otu1, otu2)) -> cor.df.ann
  
  #only use OTU pairs where both are IDd to genus
  cor.df.ann %>%
    filter(pair_haveGenusIDs == "yes") -> cor.df.ann
  
  ########################
  #enviro
  
  #identify correlations between OTU pairs that are stable
  cor.df.ann %>%
    select(source, otupair, enviro_cor, enviro_signif) %>%
    group_by(otupair) %>%
    summarize(numSignifLevels = length(unique(enviro_signif)),
              signifLevels = paste(unique(enviro_signif), collapse = "_")) -> summ.enviro
  summ.enviro$signifLevels <- recode(summ.enviro$signifLevels, negative_no = "no_negative", positive_no = "no_positive",
                                     negative_no_positive = "no_negative_positive", positive_no_negative = "no_negative_positive")
  summ.enviro %>%
    filter(signifLevels %in% c("negative","positive")) -> stable.otupairs
  # make the dataframe
  cor.df.ann %>%
    filter(source == "run1") %>%
    filter(otupair %in% stable.otupairs$otupair) %>%
    arrange(desc(enviro_cor)) -> env.frame
  
  # pull out the key OTU pairs
  lastrow <- dim(env.frame)[1]
  select.env <- env.frame[c(1:3, lastrow - 2, lastrow - 1, lastrow),]
  
  # merge with X coefs associated with each
  Xcoefs.list <- lapply(fit.list2, function(x) x$Xcoefs.df)
  
  Xcoefs.summary.list <- lapply(Xcoefs.list, function(x){
    result<- summary_OTUID_nonsignif_Xcoefs(Xcoefs.df = x, taxAndFunguild = taxAndFunguild)
    result.df<-data.frame(result)
    return(result.df)
  })
  Xcoefs.summary.df <- list_to_df(Xcoefs.summary.list)
  
  #pull out Xcoef info for each pair
  xcoef.sub.list <- list()
  for(i in 1:dim(select.env)[1]){
    
    otu1 <- select.env$otu1[i]
    otu2 <- select.env$otu2[i]
    
    Xcoefs.summary.df %>%
      filter(source == "run1") %>%
      filter(OTUId %in% c(otu1, otu2)) %>%
      select(-source) %>%
      gather(key = "key", value = "value", -c(OTUId)) %>%
      spread(key = OTUId, value = value) %>%
      mutate(enviro_cor = select.env$enviro_cor[i]) %>%
      mutate(otupair = select.env$otupair[i]) -> xcoef.sub.list[[i]]
  }
  names(xcoef.sub.list) <- select.env$otupair
  
  # negative example
  length(xcoef.sub.list)
  negExample <- xcoef.sub.list[[6]]
  taxAndFunguild %>%
    filter(OTUId %in% c(colnames(negExample)[2:3])) %>%
    select(OTUId, OTUId_ann) -> names.indx
  colnames(negExample)[2:3] <- names.indx$OTUId_ann

  modelID2 <- "allX"
  fileName <- paste0("output/boral_cooccur/chordDiagrams_enviro_negExamp", modelID2, ".csv")
  write.csv(negExample, file = fileName)

  # taxAndFunguild %>%
  #   filter(OTUId %in% c(colnames(xcoef.sub.list[[3]][2:3]))) %>%
  #   select(OTUId, OTUId_ann)
  # xcoef.sub.list[[1]]
  # xcoef.sub.list[[2]]
  # xcoef.sub.list[[3]]
  
}

compareAbund <- function(otu1, otu2, taxAndFunguild, otu.tab, covariates, seqSamples, b.levels, corval){
  
  # select otu names and abundances
  otu1.name <- taxAndFunguild[taxAndFunguild$OTUId == otu1, "OTUId_ann"]
  otu2.name <- taxAndFunguild[taxAndFunguild$OTUId == otu2, "OTUId_ann"]
  otu.tab.1 <- otu.tab[,colnames(otu.tab) == otu1]
  otu.tab.2 <- otu.tab[,colnames(otu.tab) == otu2]
  sub.otu.tab <- data.frame(seq_sampName = row.names(otu.tab), otu1 = otu.tab.1, otu2 = otu.tab.2)
  
  # add binomial data and calc mean OTUabund by code
  sub.otu.tab %>%
    gather(key = "OTUId", value = "OTUabund", -seq_sampName) %>%
    left_join(covariates) %>%
    left_join(seqSamples) %>%
    group_by(OTUId, Binomial, size) %>%
    summarize(mean_abund = round(mean(OTUabund), digits = 2)) -> abund.df
  b.levels.use <- b.levels[b.levels %in% unique(abund.df$Binomial)] # removes Olax
  abund.df$Binomial <- factor(abund.df$Binomial, levels = b.levels.use)
  
  #plot
  p <- ggplot(abund.df, aes(y = Binomial, x = OTUId)) +
    geom_tile(aes(fill = mean_abund), color = "gray") +
    geom_text(aes(label = mean_abund), size = 3) +
    facet_grid(~ size) +
    scale_fill_gradient2(name = "OTU abundance", 
                         mid = "gray",
                         na.value = "white") +
    theme_classic() + scale_x_discrete(position = "top") +
    ylab("Wood species") + xlab("") + guides(fill = FALSE) +
    ggtitle(paste0(otu1.name, " and ", otu2.name, ", residual cor =", round(corval, digits = 2)))

  return(p)
  
}

investigate_resid_chordDiagram <- function(fit.list, fit.list2, taxAndFunguild, allXs, complete_subset_list, seqSamples, zanneTree){
  
  # order the Binomial names by the phylo tree
  ggt <- ggtree(zanneTree) + geom_tiplab(size=4)
  phylo.order <- ggt$data[1:22, c("y","label")]
  phylo.order %>%
    rename("Binomial.order"="y",
           "Binomial"="label") %>%
    arrange(Binomial.order) -> phylo.order
  b.levels <- phylo.order$Binomial
  b.levels <- gsub("_"," ", b.levels)
  
  # name the otu table and covariates
  otu.tab <- complete_subset_list[['otus.trimmed']]
  covariates <- complete_subset_list[['covariates']]
  covariates %>%
    select(seq_sampName, code) -> covariates
  
  # create the correlation dataframe
  cor.list <- lapply(fit.list, function(x) x$cor.df)
  cor.df <- list_to_df(cor.list)
  cor.df.ann <- annotate_cor_withOTUInfo(cor.df = cor.df, taxAndFunguild = taxAndFunguild)
  cor.df.ann %>%
    mutate(otupair = paste0(otu1, otu2)) -> cor.df.ann
  
  #only use OTU pairs where both are IDd to genus
  cor.df.ann %>%
    filter(pair_haveGenusIDs == "yes") -> cor.df.ann
  
  ########################
  #residual
  
  #identify correlations between OTU pairs that are stable
  cor.df.ann %>%
    select(source, otupair, residual_cor, residual_signif) %>%
    group_by(otupair) %>%
    summarize(numSignifLevels = length(unique(residual_signif)),
              signifLevels = paste(unique(residual_signif), collapse = "_")) -> summ.residual
  summ.residual$signifLevels <- recode(summ.residual$signifLevels, negative_no = "no_negative", positive_no = "no_positive",
                                       negative_no_positive = "no_negative_positive", positive_no_negative = "no_negative_positive")
  summ.residual %>%
    filter(signifLevels %in% c("negative","positive")) -> stable.otupairs
  
  # make the dataframe
  cor.df.ann %>%
    filter(source == "run1") %>%
    filter(otupair %in% stable.otupairs$otupair) %>%
    arrange(desc(residual_cor)) -> res.frame
  
  # pull out the key OTU pairs
  lastrow <- dim(res.frame)[1]
  select.res <- res.frame[c(1:3, lastrow - 2, lastrow - 1, lastrow),]
  
  # pull out their OTU abund across wood species and class sizes and plot
  p.list <- list()
  for(i in 1:dim(select.res)[1]){
    p.list[[i]] <- compareAbund(otu1 = select.res$otu1[i], otu2 = select.res$otu2[i], 
                                taxAndFunguild, otu.tab, covariates, seqSamples, b.levels,
                                corval = select.res$residual_cor[i])
  }

  return(p.list)
  
  #p.list
  
  # taxAndFunguild %>%
  #   filter(OTUId == select.res$otu2[1])
  
}

make_chordDiagrams_table <- function(fit.list, taxAndFunguild, complete_subset_list, allXs = TRUE){
  
  if(allXs == TRUE){
    modelID2 <- "allX"
  }else{
    modelID2 <- "selectX"
  }
  
  # create the correlation dataframe
  cor.list <- lapply(fit.list, function(x) x$cor.df)
  cor.df <- list_to_df(cor.list)
  cor.df.ann <- annotate_cor_withOTUInfo_v(cor.df = cor.df, taxAndFunguild = taxAndFunguild)
  cor.df.ann %>%
    mutate(otupair = paste0(otu1, otu2)) -> cor.df.ann
  
  #only use OTU pairs where both are IDd to genus
  cor.df.ann %>%
    filter(pair_haveGenusIDs == "yes") -> cor.df.ann
  
  ########################
  #enviro
  
  #identify correlations between OTU pairs that are stable
  cor.df.ann %>%
    select(source, otupair, enviro_cor, enviro_signif) %>%
    group_by(otupair) %>%
    summarize(numSignifLevels = length(unique(enviro_signif)),
              signifLevels = paste(unique(enviro_signif), collapse = "_")) -> summ.enviro
  summ.enviro$signifLevels <- recode(summ.enviro$signifLevels, negative_no = "no_negative", positive_no = "no_positive",
                                     negative_no_positive = "no_negative_positive", positive_no_negative = "no_negative_positive")
  summ.enviro %>%
    filter(signifLevels %in% c("negative","positive")) -> stable.otupairs
  
  # make the plotting dataframe
  cor.df.ann %>%
    filter(source == "run1") %>%
    filter(otupair %in% stable.otupairs$otupair) %>%
    arrange(desc(enviro_cor)) %>%
    select(otu1, OTUIdann_1, phylum_1, troph_1, otu2, OTUIdann_2, phylum_2, troph_2, enviro_cor) %>%
    rename('corVal'='enviro_cor') %>%
    mutate(corType = "enviro") -> env.frame
  dim(env.frame)

  ########################
  #residual
  
  #identify correlations between OTU pairs that are stable
  cor.df.ann %>%
    select(source, otupair, residual_cor, residual_signif) %>%
    group_by(otupair) %>%
    summarize(numSignifLevels = length(unique(residual_signif)),
              signifLevels = paste(unique(residual_signif), collapse = "_")) -> summ.residual
  summ.residual$signifLevels <- recode(summ.residual$signifLevels, negative_no = "no_negative", positive_no = "no_positive",
                                       negative_no_positive = "no_negative_positive", positive_no_negative = "no_negative_positive")
  summ.residual %>%
    filter(signifLevels %in% c("negative","positive")) -> stable.otupairs
  
  # make the dataframe
  cor.df.ann %>%
    filter(source == "run1") %>%
    filter(otupair %in% stable.otupairs$otupair) %>%
    arrange(desc(residual_cor)) %>%
    select(otu1, OTUIdann_1, phylum_1, troph_1, otu2, OTUIdann_2, phylum_2, troph_2, residual_cor) %>%
    rename('corVal'='residual_cor') %>%
    mutate(corType = "residual")-> res.frame
  dim(res.frame)
  
  
  #########################
  #make table
  
  frame <- rbind(env.frame, res.frame)
  frame$phylum_1 <- recode(frame$phylum_1, Basidiomycota = "B", Ascomycota = "A")
  frame$phylum_2 <- recode(frame$phylum_2, Basidiomycota = "B", Ascomycota = "A")
  frame$troph_1 <- recode(frame$troph_1, Saprotroph = "Sapro", Pathotroph = "Patho", Symbiotroph = "Symbio")
  frame$troph_2 <- recode(frame$troph_2, Saprotroph = "Sapro", Pathotroph = "Patho", Symbiotroph = "Symbio", 
         'Pathotroph-Symbiotroph' = "Patho-Symbio")
  frame$corVal <- round(frame$corVal, digits = 2)
  
  # for each OTU, look up the proportion of reads and the number of occurences
  un.OTUId <- unique(c(frame$otu1, frame$otu2))
  comm.otu <- complete_subset_list$otus.trimmed
  result.list <- list()
  for(i in 1:length(un.OTUId)){
    curr.otu.abund <- comm.otu[,colnames(comm.otu) == un.OTUId[i]]
    numReads <- sum(curr.otu.abund)
    numSamps <- sum(curr.otu.abund > 0)
    result.list[[i]]<- data.frame(OTUId = un.OTUId[i], numReads, numSamps)
  }
  result.df <- list_to_df(result.list)
  
  colnames(result.df)[1]<-"otu1"
  frame %>%
    left_join(result.df) %>%
    rename('numReads_1'='numReads',
           'numSamps_1'='numSamps') -> tmp
  colnames(result.df)[1]<-"otu2"
  tmp %>%
    left_join(result.df) %>%
    rename('numReads_2'='numReads',
           'numSamps_2'='numSamps') -> frame.ann
  
  totalReads <- sum(colSums(comm.otu))
  frame.ann %>%
    mutate(percReads_1 = round((numReads_1 / totalReads) *100, digits = 2)) %>%
    mutate(percReads_2 = round((numReads_2 / totalReads) *100, digits = 2)) %>%
    select(OTUIdann_1, numSamps_1, percReads_1, phylum_1, troph_1, 
           OTUIdann_2, numSamps_2, percReads_2, phylum_2, troph_2, 
           corVal, corType) -> frame.final
  
  fileName <- paste0("output/boral_cooccur/chordDiagrams_", modelID2, ".csv")
  write.csv(frame.final, file = fileName)
  
}



