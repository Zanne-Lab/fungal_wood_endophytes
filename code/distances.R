# phylogenetic, functional, and community distances
#######


#---------------------------------------------------------#
# examine wood phylogenetic distance dimensions

load_zanneTreePCOs <- function(zanneTree) {
  #PCoA using cophenetic distances
  zanneTree.pco <- scores(capscale(as.dist(cophenetic(zanneTree)) ~ 1), choice=1:21, display='sites')
  return(zanneTree.pco)
}

plot_zanneTreePCO <- function(zanneTree) {
  
  require(ggtree)
  
  zanneTree.pco <- load_zanneTreePCOs(zanneTree)
  ggt <- ggtree(zanneTree) + geom_tiplab(size=4)
  graphing_df1<-data.frame(zanneTree.pco)
  p.treePCoAs<-gheatmap(ggt, graphing_df1, 
           offset = 100, width=0.7, 
           colnames_position="top",
           font.size=2, colnames_angle = 90, hjust=.2,
           low = "white", high = "blue") + ggtitle("PCoA site scores")
  
  p.treePCoAs
  ggsave("output/distances/woodPhylo_pcoaScores.pdf", width=8, height=6)
  
}

findSignif_zanneTreePCOs <- function(zanneTree, seqSamples, comm.otu){
  
  require(vegan)
  
  zanneTree.pco <- load_zanneTreePCOs(zanneTree)
  
  # link the wood species PCoAs with seqSamples
  df.xVars<-data.frame(Binomial=row.names(zanneTree.pco), zanneTree.pco)
  df.xVars %>%
    mutate(Binomial = gsub("_"," ", Binomial)) %>%
    left_join(seqSamples) %>%
    arrange(seq_sampName) -> df.xVars.samp
  dataCols<-grepl("M", colnames(df.xVars.samp))
  xVars.mat<-df.xVars.samp[,dataCols]
  row.names(xVars.mat)<-df.xVars.samp$seq_sampName
  
  #check row names in the community matrix are the same order
  row.names(xVars.mat) == row.names(comm.otu)
  o<-match(row.names(xVars.mat), row.names(comm.otu))
  comm.otu.o<-comm.otu[o,]
  sum(row.names(xVars.mat) != row.names(comm.otu.o)) # this needs to be 0
  
  # use dbRDA to evaluate which PCoA axes affect community composition
  # where matrix includes site and size
  cap.env <- capscale(comm.otu.o ~ ., data=xVars.mat, distance='bray')
  # set up the null cases with no predictors
  mod0.env <- capscale(comm.otu.o ~ 1, data=xVars.mat, distance='bray')
  # model selection
  step.env <- ordistep(mod0.env, scope=formula(cap.env))
  xs<-unlist(lapply(strsplit(row.names(step.env$anova), "+ "), function(x){x[2]}))
  dbrda.anova <- data.frame(terms=xs, step.env$anova)

  #all are significant in the model
  #some of the early splits may be slightly more informative, but these are not dramatically so
  #pretty much the same as putting species in the model
  
  write.csv(dbrda.anova, file="output/distances/woodPhylo_dbrdaANOVA.csv", row.names=FALSE)
  
}


#---------------------------------------------------------#
# make wood trait ordination

plot_woodTrait_pca <- function(traits.code, seqSamples){
  
  require(ggplot2)
  
  # calculate wood functional trait distance in multivariate space 
  
  # identify rows with no missing values
  x <- complete.cases(traits.code[,-(1:3)]) 
  traits.mean1<-traits.code[x,-(1:3)]
  traits.samps<-traits.code[x,1:3]
  
  # match wood species to families
  famIndx<-unique(seqSamples[,c("species","Family")])
  traits.samps %>%
    left_join(famIndx) -> traits.samps
  # Olax stricta is missing because it doesn't have a waterperc value and it is only represented in small stem samples 
  
  #log-transform and scale, do PCA and take the first 2 axes
  traits.scaled <- apply(log(traits.mean1+10), 2, scale)
  pc <- prcomp(traits.scaled)
  df <- data.frame(pc$x[,1:2], traits.samps)
  df$size1<-factor(df$size)
  
  #put together the loadings dataframe
  datapc <- data.frame(varnames=rownames(pc$rotation), pc$rotation[,1:2])
  mult <- min(
    (max(df[,"PC1"]) - min(df[,"PC2"])/(max(datapc[,"PC2"])-min(datapc[,"PC2"]))),
    (max(df[,"PC2"]) - min(df[,"PC1"])/(max(datapc[,"PC1"])-min(datapc[,"PC1"])))
  )
  datapc <- transform(datapc,
                      v1 = .7 * mult * datapc$PC1,
                      v2 = .7 * mult * datapc$PC2)
  #plot
  mytheme <- make_ggplot_theme()
  plot <- ggplot(data=datapc, aes(x=v1, y=v2, label=varnames)) + 
    geom_text(size=3, vjust=1, color=1) + 
    geom_segment(data=datapc, aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color=1) + 
    coord_fixed() + 
    geom_point(mapping=aes(x=PC1, y=PC2, color=Family, shape=size1), data=df, inherit.aes = FALSE) +
    mytheme + xlab("PC1 (24%)") + ylab("PC2 (18%)") +
    scale_shape_manual(name="Size class", values=c(16,17)) + 
    scale_color_viridis(name="Wood family", discrete=TRUE)
  
  ggsave("output/distances/PCAofWoodTraits.pdf", width=5, height=4)
  
}


#---------------------------------------------------------#
# calculate distances

# extract_uniquePairDists() is specified in load_data.R

calc_woodPhyloDist <- function(zanneTree){
  
  # calc wood phylogenetic distance
  phyloDist <- cophenetic(zanneTree)
  mat.phyloDist <- as.matrix(phyloDist)
  
  #make it long
  phyloDist.l <- extract_uniquePairDists(dist.mat=mat.phyloDist) 
  colnames(phyloDist.l)<-c("samp1_sp","samp2_sp","woodPhyloDist")
  
  #check the number of unique species
  allSp <- unique(c(as.character(phyloDist.l$samp1_sp), 
                    as.character(phyloDist.l$samp2_sp)))
  #length(allSp) #should be 22
  #add species vs itself
  sameSp<-data.frame(samp1_sp = allSp, 
                     samp2_sp = allSp, 
                     woodPhyloDist=rep(0,length(allSp))
                     )
  phyloDist.l<-rbind(phyloDist.l,sameSp)
  
  phyloDist <- phyloDist.l
  return(phyloDist)
}

calc_woodTraitDist <- function(traits.code){
  
  # identify rows with no missing values
  traits.code %>%
    filter(complete.cases(waterperc, density, barkthick, P, K, Ca, Mn, Fe, Zn, N, C)) -> trait.code.c
  trait.code.c %>%
    select(waterperc, density, barkthick, P, K, Ca, Mn, Fe, Zn, N, C) -> trait.code.mat
  # this step filtered out Olax stricta because it doesn't have a waterperc value and it is only represented in small stem samples 
  
  #log-transform and scale, do PCA and take the first 3 axis
  traits.scaled <- apply(log(trait.code.mat+10), 2, scale)
  pc <- princomp(traits.scaled)  # stats package
  pc.scores <- pc$scores[, 1:3]  # the first 3 axes
  
  # make a unique identifier for each row
  row.names(pc.scores) <- trait.code.c$code
  pc.dist.mat <- dist(pc.scores, method = "euclidean", diag=TRUE, upper=TRUE) #calc euclidean distance
  mat.traitDist<-as.matrix(pc.dist.mat)
  
  #make it long
  traitDist.l <- extract_uniquePairDists(dist.mat=mat.traitDist) #make it long
  colnames(traitDist.l)<-c("code1","code2","woodTraitDist")
  
  #add code vs itself
  allCodes<-unique(c(as.character(traitDist.l$code1), 
                     as.character(traitDist.l$code2))) # 32 codes
  sameCode<-data.frame(code1 = allCodes,
                       code2 = allCodes,
                       woodTraitDist = 0)
  
  traitDist.l<-rbind(traitDist.l, sameCode)

  traitDist <- traitDist.l
  return(traitDist)
}

calc_commDist <- function(comm.otu){
  
  require(vegan)
  
  # calc fungal community dissimilarity distance
  fungDist <-vegdist(decostand(comm.otu, 'hellinger'),'bray') #calculate bray-curtis distance
  mat.fungDist <- as.matrix(fungDist)
  
  #make it long
  fungDist.l <- extract_uniquePairDists(dist.mat=mat.fungDist) #make it long
  colnames(fungDist.l)<-c("samp1","samp2","commDist")
  
  commDist <- fungDist.l
  return(commDist)
}

calc_pairdist <- function(zanneTree, comm.otu, traits.code, seqSamples){
  
  #calculate distances
  phyloDist <- calc_woodPhyloDist(zanneTree)
  commDist <- calc_commDist(comm.otu)
  traitDist <- calc_woodTraitDist(traits.code)
  
  #join phyloDist (pairs of Binomial) and commDist (pairs of seq_sampName)
  #annotate commDist samples with Binomial data
  seqSamples %>%
    select(seq_sampName, Binomial) %>%
    mutate(Binomial = gsub(" ", "_", Binomial)) -> binom.indx
  #merge col1
  binom.indx %>%
    rename(samp1=seq_sampName,
           samp1_sp=Binomial) -> binom.indx
  commDist %>%
    left_join(binom.indx) -> commDist.1
  #merge col2
  binom.indx %>%
    rename(samp2=samp1,
           samp2_sp=samp1_sp) -> binom.indx
  commDist.1 %>%
    left_join(binom.indx) -> commDist.2
  #find unique Binomial pairs in commDist.2 and phyloDist
  commDist.2 %>%
    mutate(spPair_rev = paste(samp2_sp, samp1_sp)) -> commDist.2
  phyloDist %>%
    mutate(spPair = paste(samp1_sp, samp2_sp)) -> phyloDist
  commDist.2 %>%
    left_join(phyloDist) %>%
    select(-spPair) %>%
    left_join(phyloDist, by=c("spPair_rev"="spPair")) %>%
    mutate(woodPhyloDist = ifelse(is.na(woodPhyloDist.x), woodPhyloDist.y, woodPhyloDist.x)) %>%
    select(samp1, samp2, commDist, woodPhyloDist) -> commDist.woodDist
  
  #join commDist.woodDist (pairs of seq_sampName) and traitDist (pairs of codes)
  #annotate commDist.woodDist with codes
  seqSamples %>%
    select(seq_sampName, code) -> code.indx
  #merge col1
  code.indx %>%
    rename(samp1=seq_sampName,
           code1=code) -> code.indx
  commDist.woodDist %>%
    left_join(code.indx) -> commDist.woodDist.1
  #merge col2
  colnames(code.indx)
  code.indx %>%
    rename(samp2=samp1,
           code2=code1) -> code.indx
  commDist.woodDist.1 %>%
    left_join(code.indx) -> commDist.woodDist.2
  #find unique code pairs in commDist.woodDist2 and traitDist
  commDist.woodDist.2 %>%
    mutate(codePair_rev = paste(code2, code1)) -> commDist.woodDist.2
  traitDist %>%
    mutate(codePair = paste(code1, code2)) -> traitDist
  commDist.woodDist.2 %>%
    left_join(traitDist) %>%
    select(-codePair) %>%
    left_join(traitDist, by=c("codePair_rev"="codePair")) %>%
    mutate(woodTraitDist = ifelse(is.na(woodTraitDist.x), woodTraitDist.y, woodTraitDist.x)) %>%
    select(samp1, samp2, commDist, woodPhyloDist, woodTraitDist) -> commDist.woodDist.traitDist
 
  #annotate with site, size, binomial
  fam.indx<-unique(seqSamples[,c("code","Family")])
  seqSamples %>%
    left_join(fam.indx) %>%
    select(seq_sampName, species, size, Binomial, site, Family) %>%
    separate(Binomial, into=c("genus","speciesName")) %>%
    select(-speciesName) %>%
    rename(samp1=seq_sampName) -> seqIndx
  commDist.woodDist.traitDist %>%
    left_join(seqIndx) %>%
    rename(samp1_species=species,
           samp1_size=size,
           samp1_genus=genus,
           samp1_site=site,
           samp1_family=Family) -> commDist.woodDist.traitDist.1
  seqIndx %>%
    rename(samp2=samp1) -> seqIndx
  commDist.woodDist.traitDist.1 %>%
    left_join(seqIndx) %>%
    rename(samp2_species=species,
           samp2_size=size,
           samp2_genus=genus,
           samp2_site=site,
           samp2_family=Family) -> commDist.woodDist.traitDist.2

  # prune and annotate particular sample pairs
  commDist.woodDist.traitDist.2 %>%
    filter(samp1_size == samp2_size) %>% # remove pairwise comparision between large and small stem size samples
    mutate(sameSite = ifelse(samp1_site == samp2_site, "yes", "no")) %>%
    mutate(sameSpecies = ifelse(samp1_species == samp2_species, "yes", "no")) %>%
    mutate(sameGenus = ifelse(samp1_genus == samp2_genus, "yes","no")) %>%
    mutate(sameFamily = ifelse(samp1_family == samp2_family, "yes", "no")) -> pairDist

  #add category names
  pairDist$woodPhyloCat <- "betweenFamilies" #categorize level of phylogenetic similarity
  pairDist[pairDist$sameSpecies == "yes", "woodPhyloCat"] <-"sameSpecies"
  pairDist[pairDist$sameSpecies == "no" & pairDist$sameGenus == "yes", "woodPhyloCat"] <- "sameGenus"
  pairDist[pairDist$sameSpecies == "no" &
             pairDist$sameGenus == "no" &
             pairDist$sameFamily == "yes", "woodPhyloCat"] <- "sameFamily"
  pairDist$woodPhyloCat <-
    factor(
      pairDist$woodPhyloCat,
      levels = c("sameSpecies", "sameGenus", "sameFamily", "betweenFamilies")
    )
  
  # prune cols
  pairDist %>%
    select(samp1, samp2, commDist, woodPhyloDist, woodTraitDist, samp1_size, samp1_species, samp2_species, sameSite, woodPhyloCat) %>%
    rename(size = samp1_size) -> pairDist #removed pairs between sizes, so this is ok
  
  return(pairDist)
}


#---------------------------------------------------------#
# examine wood phylogenetic vs fungal community distances

plot_woodDistfungDist_bothContinuous <- function(pairDist){
  
  require(ggplot2)
  require(viridis)
  require(gridExtra)
  
  mytheme <- make_ggplot_theme()
  
  p1.phylo <- ggplot(pairDist, aes(x = woodPhyloDist, y = commDist, color=size)) +
    geom_point(alpha=0.25, pch=16) + geom_smooth(formula=y~poly(x,2), method="lm") +
    xlab("Wood phylogenetic distance") +
    ylab("Fungal community distance") + mytheme + guides(color=FALSE) +
    scale_color_viridis(discrete = TRUE, option ="D", end=0.8)
  
  p2.trait<-ggplot(pairDist, aes(x=as.numeric(woodTraitDist), y=commDist, color=size)) + 
    geom_point(alpha=0.25, pch=16) + 
    xlab("Wood trait distance") + 
    ylab("Fungal community distance") + mytheme + guides(color=FALSE) +
    scale_color_viridis(discrete = TRUE, option ="D", end=0.8) + 
    geom_smooth(formula=y~poly(x,2), method="lm")
  
  #remove 0,0 points
  temp<-pairDist[pairDist$woodPhyloDist!=0 & !is.na(pairDist$woodPhyloDist) &
                   pairDist$woodTraitDist!=0 & !is.na(pairDist$woodTraitDist),]

  #use only unique wood species x size (ie not the reps)
  temp$uniqTrtPair<-paste(temp$samp1_sp, temp$samp2_sp, temp$size, sep="_")
  temp.df<-unique(temp[,c("uniqTrtPair","woodPhyloDist","woodTraitDist","size")])

  p3.phyloTrait<-ggplot(temp.df, aes(x=woodPhyloDist, 
                                  y=as.numeric(woodTraitDist), color=size)) + 
    geom_point(alpha=0.5) + 
    xlab("Wood phylogenetic distance") + 
    ylab("Wood trait distance") + mytheme + 
    guides(color=FALSE) +
    scale_color_viridis(name="Size class",discrete = TRUE, option ="D", end=0.8)
  
  #extract legend
  p<-ggplot(pairDist, aes(x=woodPhyloCat, y=as.numeric(woodTraitDist), fill=size)) + 
    geom_violin(alpha=0.9) + mytheme + scale_fill_viridis(name="Size class", discrete = TRUE, option ="D", end=0.8)
  legend <- g_legend(p) 
  
  pdf(file="output/distances/fig_woodPhylo_woodTrait_fungDist.pdf", height=6, width=6, colormodel = "cmyk")
  grid.arrange(p2.trait+ggtitle("a"), p1.phylo+ggtitle("b"),
               p3.phyloTrait+ggtitle("c"), legend, ncol=2)
  dev.off()
  
}

write_woodDistfungDist_summary<-function(pairDist){
  
  #phylo
  mod<-lm(commDist ~ woodPhyloDist + size + woodPhyloDist:size, data=pairDist) 
  mod2.phylo<-lm(commDist ~ woodPhyloDist + size + woodPhyloDist:size + I(woodPhyloDist^2), data=pairDist)
  #anova(mod, mod2.phylo)# squared term works better
  phylo.sum<-summary(mod2.phylo)$coefficients
  
  #functional traits
  mod<-lm(commDist ~ woodTraitDist + size + woodTraitDist:size, data=pairDist) 
  mod2.funct<-lm(commDist ~ woodTraitDist + size + woodTraitDist:size + I(woodTraitDist^2), data=pairDist)
  #anova(mod, mod2.funct)# squared term works better
  funct.sum<-summary(mod2.funct)$coefficients

  #make the table
  tab1<-rbind(funct.sum, phylo.sum)
  
  write.csv(tab1, file="output/distances/summary_woodPhylo_woodTrait_fungDist.csv", row.names=TRUE)
}




