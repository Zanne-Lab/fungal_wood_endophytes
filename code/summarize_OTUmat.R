# summarize the OTU matrix
#######

calc_otu_summStats<-function(comm.otu, taxAndFunguild){
  
  require(vegan)
  
  #total number of OTUs
  totalOTUs<-dim(comm.otu)[2]
  
  #mean sample richness
  richness<-apply(comm.otu, MARGIN = 1,function(x) sum(x>0))
  meanRichness<-mean(richness)
  seRichness<-sd(richness)/sqrt(length(richness))
  
  # #mean sample evenness
  # H <- diversity(comm.otu) 
  # S <- specnumber(comm.otu) ## rowSums(BCI > 0) does the same...
  # J <- H/log(S) #Pielou's evenness (J)
  # meanJ<-mean(J)
  # seJ<-sd(J)/sqrt(length(J))
  
  #mean number of reads per sample
  meanReads<-mean(rowSums(comm.otu))
  seReads<-sd(rowSums(comm.otu))/sqrt(length(rowSums(comm.otu)))
  
  #taxon ID coverage
  taxAndFunguild %>%
    select(OTUId, genus, phylum, Trophic.Mode) %>%
    filter(OTUId %in% colnames(comm.otu)) -> tax.indx
  tax.indx %>%
    filter(genus == "unclassified") -> genus.unclass
  num.genus.unclass <- dim(genus.unclass)[1]
  tax.indx %>%
    filter(genus != "unclassified") %>%
    filter(Trophic.Mode == "unclassified") -> guild.unclass
  num.guild.unclass <- dim(guild.unclass)[1]
  
  #what % of OTUs are present in abundances > 1% of total sequence reads?
  totalSeq <- sum(rowSums(comm.otu))
  oneperc.totalSeq <- totalSeq /100
  numOTUsgreater <- sum(colSums(comm.otu) > oneperc.totalSeq)
  numOTUsgreater
  
  #what % of OTUs are observed in fewer than X samples?
  comm.otu.pa <- (comm.otu > 1) *1
  numOTUsfewer<- sum(colSums(comm.otu.pa) < 10)
  (numOTUsfewer / totalOTUs) *100
  
    summaryStats<-data.frame(label=c("totalOTUs","meanRichness","seRichness","meanReads","seReads", "num.genus.unclass","num.guild.unclass"),
                           value=c(totalOTUs, meanRichness, seRichness, meanReads, seReads, num.genus.unclass, num.guild.unclass))
  summaryStats$value <- round(summaryStats$value, digits=4)

  write.csv(summaryStats, file="output/otuSummary/otuSummary.csv")
}

plot_sampleEffortCurves<-function(comm.otu){
  
  require(vegan)
  
  pdf(file="output/otuSummary/sampleEffortCurves.pdf", width=5, height=5)
  rarecurve(comm.otu, step=100, xlab="Number of reads per sample", ylab="Cumulative number of OTUs", label=FALSE)
  dev.off()
  
}

plot_fungal_richness_plant_phylo<-function(zanneTree, comm.otu, seqSamples){
  
  require(ggtree)
  require(viridis)
  
  #plot wood phylo
  ggt <- ggtree(zanneTree) + geom_tiplab(size=4)
  
  #calc OTU richness in each sample, average over code
  richness <- apply(comm.otu, MARGIN = 1, function(x) sum(x>0)) 
  richness_df <- data.frame(seq_sampName=names(richness), richness=richness)
  richness_df %>%
    left_join(seqSamples) %>%
    group_by(Binomial, size) %>%
    summarize(mean_otu_rich = mean(richness)) %>%
    spread(key=size, value=mean_otu_rich) -> species_mean_df
  
  #pretty graphing df
  graphing_df<-data.frame(Small=species_mean_df$small, Large=species_mean_df$large)
  row.names(graphing_df)<-sub(" ","_",species_mean_df$Binomial)
  gheatmap(ggt, graphing_df, 
           offset = 100, width=0.4, 
           colnames_position="top",
           font.size=4) +
    viridis::scale_fill_viridis(option="D",na.value = "white")
  
  ggsave("output/otuSummary/fungal_richness_fig.pdf", width=7, height=6)
}

calc_subsetOTU_richness <- function(otuids, comm.otu){
  
  comm.otu.subset <- comm.otu[,colnames(comm.otu) %in% otuids]
  comm.otu.subset.pa <- (comm.otu.subset > 0) * 1 # turn table into presence/absence
  df <- data.frame(seq_sampName = names(rowSums(comm.otu.subset.pa)), 
             otusub.rich = rowSums(comm.otu.subset.pa))
  return(df)

  }

plot_richness <- function(comm.otu, seqSamples, taxAndFunguild, zanneTree){
  
  require(ggtree)
  
  # oomycete richness
  # taxAndFunguild %>%
  #   filter(kingdom == "Protist") %>%
  #   select(OTUId) -> sub.otuids
  # oo.df <- calc_subsetOTU_richness(otuids = sub.otuids$OTUId, comm.otu)
  # colnames(oo.df)[2] <- "Oomycetes"

  # fungal richness
  taxAndFunguild %>%
    filter(kingdom == "Fungi") %>%
    select(OTUId) -> sub.otuids
  fung.df <- calc_subsetOTU_richness(otuids = sub.otuids$OTUId, comm.otu)
  colnames(fung.df)[2] <- "Fungi"
  
  # combine and plot
  fung.df %>%
    gather(key = "type", value ="richness", -1) %>%
    left_join(seqSamples) %>%
    group_by(Binomial, size, type) %>%
    summarize(mean_rich = round(mean(richness), digits = 0)) -> tmp
  
  # order the size levels
  tmp$size<-recode(tmp$size, `small`="Small", `large`="Large")
  tmp$size<-factor(tmp$size, levels = c("Small","Large"))
  
  # order the Binomial names by the phylo tree
  ggt <- ggtree(zanneTree) + geom_tiplab(size=4)
  allSp <- unique(tmp$Binomial)
  phylo.order <- ggt$data[1:length(allSp), c("y","label")]
  phylo.order %>%
    rename("Binomial.order"="y",
           "Binomial"="label") %>%
    arrange(Binomial.order) -> phylo.order
  b.levels <- phylo.order$Binomial
  b.levels <- gsub("_"," ", b.levels)
  tmp$Binomial <- factor(tmp$Binomial, levels = b.levels)
  
  #plot
  p <- ggplot(tmp, aes(y = Binomial, x = size)) +
    geom_tile(aes(fill = mean_rich), color = "gray") +
    geom_text(aes(label = mean_rich), size = 3) +
    scale_fill_gradient2(name = "OTU richness", 
                         mid = "gray",
                         na.value = "white") +
    theme_classic() + scale_x_discrete(position = "top") +
    ylab("Wood species") + xlab("Stem size") + guides(fill = FALSE)
  p
  ggsave(filename = "output/otuSummary/richness.pdf", width=3, height=6)
  
  }

table_commonRare_otus <- function(comm.otu, seqSamples, taxAndFunguild){
  
  comm.otu.pa <- (comm.otu > 0) * 1 # turn table into presence/absence
  rare.df <- data.frame(OTUId = names(colSums(comm.otu.pa)), 
             n.samps.pres = colSums(comm.otu.pa))
  rare.df %>%
    left_join(taxAndFunguild) %>%
    mutate(basidio = grepl("Basidiomycota", phylum)) %>%
    mutate(sapro = grepl("Saprotroph", Trophic.Mode)) %>%
    mutate(patho = grepl("Pathotroph", Trophic.Mode)) %>%
    select(OTUId, OTUId_ann, n.samps.pres, basidio, sapro, patho) %>%
    arrange(desc(n.samps.pres), basidio, sapro, patho) -> plot.df
  plot.df$rare.rank <- 1:dim(plot.df)[1]
  
  # all
  p <- ggplot(plot.df, aes(x = rare.rank, y = n.samps.pres)) +
    geom_bar(stat = "identity") +
    xlab("OTU") + ylab("Number of samples")
  
  # most common
  plot.df %>%
    filter(n.samps.pres > 20) %>%
    mutate(sapro_patho = paste0(sapro, patho)) -> plot.df.common
  plot.df.common$sapro_patho<- recode_factor(plot.df.common$sapro_patho, 
                FALSEFALSE = "neither",
                TRUEFALSE = "saprotroph",
                FALSETRUE = "pathotroph",
                TRUETRUE = "sapro-pathotroph")
  
  color.vec <- c("gray","green","red","brown")
  names(color.vec) <- c("neither","saprotroph","pathotroph","sapro-pathotroph")
  
  p.common <- ggplot(plot.df.common, aes(x = rare.rank, y = n.samps.pres, fill = sapro_patho)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(name = "Trophic mode", values = color.vec) +
    xlab("OTU") + ylab("Number of samples")
  
  #plot.df.common[1:50,"OTUId_ann"]
  
  pdf(file = "output/otuSummary/otu_rarity.pdf", width = 6, height = 6)
  grid.arrange(p, p.common)
  dev.off()
}