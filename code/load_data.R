# load data and general helper functions
#######


#---------------------------------------------------------#
# general helper functions

list_to_df <- function(mylist){
  
  # make a vector of row ids that correspond to the list names
  rowid.indx <- lapply(mylist, function(x) dim(x)[1])
  sourceVec.list <- list()
  for(i in 1:length(rowid.indx)){
    sourceName <- names(rowid.indx)[i]
    numRows <- rowid.indx[[i]]
    sourceVec.list[[i]] <- rep(sourceName, numRows)
  }
  rowVec <- unlist(sourceVec.list)
  
  # combine into df
  df <- data.frame(do.call(rbind, mylist))
  df$source <- rowVec
  
  return(df)
}

extract_uniquePairDists<-function(dist.mat){
  
  x<-dist.mat
  rowCol <- expand.grid(rownames(x), colnames(x))
  labs <- rowCol[as.vector(upper.tri(x,diag=F)),]
  df <- cbind(labs, x[upper.tri(x,diag=F)])
  colnames(df) <- c("sp1","sp2","dist")
  
  return(df)
}


#---------------------------------------------------------#
# load microbial community data

# add_oomycetes<-function(fung.otu){
#   
#   # read in OTU table (uclust output) and convert to matrix (rows=samples, columns=OTUs)
#   data.otu <- read.csv('data/sequencing_T0/OTUtable_oomycetes_20171020.csv', row.names=1)
#   data.df<-data.frame(seqSamp=row.names(data.otu), data.otu)
#   
#   #make mat.otu of fungal taxa a dataframe
#   fung.df<-data.frame(seqSamp=row.names(fung.otu), fung.otu)
#   
#   #merge by seqSamp
#   comm.df<-left_join(fung.df, data.df)
#   
#   #make NAs into 0s
#   comm.df[is.na(comm.df)]<-0
#   
#   #make dataframe into a matrix again
#   row.names(comm.df)<-comm.df$seqSamp
#   comm.mat<-comm.df[,-1]
#   comm.mat<-as.matrix(comm.mat)
#   
#   return(comm.mat)
# }

load_matotu<-function(){
  
  # read in OTU table (uclust output) and convert to matrix (rows=samples, columns=OTUs)
  data.otu <- read.csv('data/sequencing_T0/DP16_OTUtable.csv', stringsAsFactors = FALSE)
  mat.otu <- as.matrix(data.otu[, 2:ncol(data.otu)]); rownames(mat.otu) <- data.otu[, 1]
  # remove extra blank sample
  mat.otu <- mat.otu[-which(rownames(mat.otu) == 'X..blank'), ]
  mat.otu <- mat.otu[, colSums(mat.otu) > 0]
  
  sum(colSums(mat.otu)==0) # if 0, then there are no empty columns
  
  # read in dataframe that contains sample information (also used to create meta and xrf)
  data <- read.csv('data/sequencing_T0/NextGenSeqencing_2016_Sample_MasterSpreadsheet.csv', stringsAsFactors=F)
  
  # re-label rownames in 'mat.otu' with sample codes
  rownames(mat.otu) <- gsub('.', '-', rownames(mat.otu), fixed=T)
  rownames(mat.otu)[match(data$NextGenID, rownames(mat.otu))] <- data$SampleCode
  
  
  # extract 'blank' and 'mock' samples from 'mat.otu', delete from 'data'
  blank <- mat.otu[grep('blank', rownames(mat.otu)), ]
  mock <- mat.otu['mock', ]
  mat.otu <- mat.otu[-c(grep('blank', rownames(mat.otu)), grep('mock', rownames(mat.otu))), ]
  
  # otus, taxa in mock (select cut-off of >=9 reads in a sample)
  tax <-read.delim('data/sequencing_T0/DP16_tax.txt', stringsAsFactors = F)
  mock <- data.frame(reads=sort(mock[mock > 0]))
  mock <- cbind(mock, tax[match(rownames(mock), tax$qseqid), 'species'])
  #mock
  mat.otu[mat.otu < 9] <- 0
  
  # otus, taxa in blank
  blank <- data.frame(reads=sort(blank[blank > 0]))
  blank <- cbind(blank, tax[match(rownames(blank), tax$qseqid), 'species'])
  #mat.otu[,'ITSall_OTUd_3713'] # the most abundant OTU in the blank does not show up in any of the samples
  
  # re-order rows in 'mat.otu' to match rows in 'data' after deleting 'blank' from 'data'
  data <- data[!data$SampleCode == 'blank', ]
  mat.otu <- mat.otu[match(data$SampleCode, rownames(mat.otu)), ]
  all(rownames(mat.otu) == data$SampleCode)  # is TRUE
  
  # add oomycetes
  #mat.otu <- add_oomycetes(mat.otu)
  
  #trim OTUs that do not show up in ANY of the samples
  x <- apply(mat.otu > 0, 2, sum)
  mat.otu <- mat.otu[, x >= 1]
  
  return(mat.otu)
  
}

load_TaxAndFunguild <- function(comm.otu.tmp){
  
  # load fungal OTU info
  funguild <-read.delim('data/sequencing_T0/DP16_funguild.txt', stringsAsFactors = F)
  tax <-read.delim('data/sequencing_T0/DP16_tax.txt', stringsAsFactors = F)
  
  # merge the dataframes by OTUId
  colnames(tax)[1] <- "OTUId" #make this match the column name in funguild
  taxAndFunguild <- left_join(tax, funguild)
  
  # delete OTUs from taxAndFunguild if not found in comm.otu (only found in blanks, mock, or very infrequently)
  # also delete OTUs with suspect coverage (probably nonfungal)
  # hist(taxAndFunguild$coverage)
  taxAndFunguild <- filter(taxAndFunguild, OTUId %in% colnames(comm.otu.tmp) & coverage > 0.9)
  # hist(taxAndFunguild$coverage)
  
  # delete OTUs from comm.otu not found in taxAndFunguild (probably plant DNA)
  # hist(rowSums(comm.otu))
  comm.otu.tmp <- comm.otu.tmp[, colnames(comm.otu.tmp) %in% taxAndFunguild$OTUId]
  # hist(rowSums(comm.otu))
  
  # # create a kingdom column
  # taxAndFunguild$kingdom<-NA
  # taxAndFunguild[grepl("Fungi", taxAndFunguild$taxonomy),"kingdom"]<-"Fungi"
  
  # add oomycete OTUs as rows
  #ooOTUs<-colnames(comm.otu)[!colnames(comm.otu) %in% taxAndFunguild$OTUId]
  #oo.df<-data.frame(OTUId=ooOTUs, kingdom="Protist")
  #taxAndFunguild<-bind_rows(taxAndFunguild, oo.df)
  
  # reorder taxAndFunguild to make OTU table
  o<-match(colnames(comm.otu.tmp), taxAndFunguild$OTUId)
  o.taxAndFunguild<-taxAndFunguild[o,]
  sum(o.taxAndFunguild$OTUId != colnames(comm.otu.tmp)) #this need to be 0
  
  # only use FUNGuild info with confidence ranking of Probable or Highly Probable
  o.taxAndFunguild[!o.taxAndFunguild$Confidence.Ranking %in% c("Probable","Highly Probable"),c("Trophic.Mode","Guild")]<-"unclassified"
  
  # select cols 
  o.taxAndFunguild %>%
    select(OTUId, taxonomy, kingdom, phylum, family, genus, species, 
           Trophic.Mode, Guild) -> o.taxAndFunguild
  
  #clean Trophic.Mode
  #unique(o.taxAndFunguild$Trophic.Mode)
  
  #clean Guild
  #unique(o.taxAndFunguild$Guild)
  o.taxAndFunguild[o.taxAndFunguild$Guild=="NULL","Guild"]<-"unclassified"
  
  #clean oomycetes
  #o.taxAndFunguild[o.taxAndFunguild$kingdom=="Protist", c("taxonomy","phylum","family","genus")]<-"unclassified"
  #o.taxAndFunguild[o.taxAndFunguild$kingdom=="Protist", c("species")]<-"unclassified_Protist"
  
  #clean taxa
  o.taxAndFunguild[is.na(o.taxAndFunguild$phylum), 'phylum'] <- 'unclassified'
  o.taxAndFunguild[is.na(o.taxAndFunguild$family), 'family'] <- 'unclassified'
  o.taxAndFunguild[is.na(o.taxAndFunguild$genus), 'genus'] <- 'unclassified'
  o.taxAndFunguild[is.na(o.taxAndFunguild$species), 'species'] <- 'unclassified'
  
  #clean species
  #species column should have [genus]_sp if the genus is known
  criteria <- o.taxAndFunguild$genus != "unclassified" & o.taxAndFunguild$species == "unclassified"
  o.taxAndFunguild %>%
    mutate(species_fake = paste(genus, "sp", sep = "_")) %>%
    mutate(species_new = ifelse(genus != "unclassified" & species == "unclassified", 
                                species_fake, species)) %>%
    select(-species_fake) %>%
    select(-species) %>%
    rename('species'='species_new') -> o.taxAndFunguild
  
  #fix weird characters in 'Montagnula_aloës'
  o.taxAndFunguild[grep('Montagnula', o.taxAndFunguild$species), "species"] <- "Montagnula_aloes"
  #add numbers to repeated names in species to indicate that they are different OTUs
  o.taxAndFunguild %>%
    filter(grepl("_sp", species)) %>%
    group_by(species) %>%
    summarize(n = length(species)) %>%
    filter(n > 1) -> spfake_indx
  for(i in 1:dim(spfake_indx)[1]){
    curr.sp <- as.character(spfake_indx[i,"species"])
    sp.vec <- o.taxAndFunguild[o.taxAndFunguild$species == curr.sp, "species"]
    new.sp.vec <- paste(sp.vec, 1:length(sp.vec), sep="")
    o.taxAndFunguild[o.taxAndFunguild$species == curr.sp, "species"] <- new.sp.vec
  }
  #simplify the OTUId and create an annotated OTUId column using species
  o.taxAndFunguild %>%
    separate(OTUId, into=c("drop","drop1","OTUId_num"), remove = FALSE) %>%
    select(-drop) %>% select(-drop1) %>%
    mutate(OTUId_simp = paste("OTU", OTUId_num, sep="_")) %>%
    select(-OTUId_num) %>%
    mutate(OTUId_ann = ifelse(species == "unclassified",
                              OTUId_simp, species)) %>%
    select(-OTUId_simp) -> o.taxAndFunguild
  
  return(o.taxAndFunguild)
}

clean_comm<-function(comm.otu.tmp, taxAndFunguild){
  # delete OTUs from comm.otu not found in taxAndFunguild (probably plant DNA)
  # hist(rowSums(comm.otu))
  comm.otu <- comm.otu.tmp[, colnames(comm.otu.tmp) %in% taxAndFunguild$OTUId]
  # hist(rowSums(comm.otu))
  
  return(comm.otu)
}


load_stemSamples<-function(){
  
  require(dplyr)
  require(tidyr)
  
  deployment <- read.csv("data/deployment.csv")
  deployment<-rename(deployment, "code"="species") #Code column
  
  #summarize by stem
  deployment %>%
    group_by(code, Stem) %>%
    summarize(num.unique=length(unique(unique))) %>%
    mutate(codeStem=paste(code, Stem, sep="")) %>%
    mutate(species=tolower(code)) %>%
    mutate(size=ifelse(code == tolower(code), 'small','large')) -> deploy.new
  
  #add species info
  species <- read.csv("data/species.csv", stringsAsFactors = FALSE)
  species[species$Binomial == "Ricinocarpus pinifolius", "Binomial"] <- "Ricinocarpos pinifolius" # fix a misspelled name
  species %>%
    mutate(species=tolower(Code)) %>%
    rename(site=Collection.site) %>%
    select(species, Family, Binomial, site) -> species.new
  deploy.new %>%
    left_join(species.new) -> stemSamples
  
  return(stemSamples)  
  
}

load_seqSamples<-function(mat.otu){
  
  stemSamples <- load_stemSamples()
  
  #identify sequence sampleIDs
  seq_indx<-data.frame(seq_sampName=row.names(mat.otu))
  
  #merge by codeStem
  stem.indx<-stemSamples[,c("codeStem","code","Stem")]
  left_join(seq_indx, stem.indx, by=c("seq_sampName"="codeStem")) %>%
    mutate(codeStem = ifelse(!is.na(Stem), paste(code, Stem, sep=""), NA)) %>%
    select(seq_sampName, codeStem) -> seqSamples.tmp
  
  #add back in code-level information for seq_samples that have been pooled by code
  code.indx <- unique(stemSamples[,c("code","species","size")])
  seqSamples.tmp %>%
    separate(seq_sampName, into=c("code","extra"), 4, remove=FALSE) %>%
    left_join(code.indx) %>%
    select(-extra) -> seqSamples
  
  #add binomial, Family, and site
  sp.indx <- unique(stemSamples[,c("species","Binomial", "Family", "site")])
  seqSamples %>%
    left_join(sp.indx) -> seqSamples
  
  return(seqSamples)
}


#---------------------------------------------------------#
# load wood phylogentic data

fix_problem_species<-function(tree, prob_species, dontchoose = wood_names){
  
  require(taxonlookup)
  
  for(i in 1:length(prob_species)){
    genus<-taxonlookup:::split_genus(prob_species[i])
    replace<-sample(tree$tip.label[grepl(genus,tree$tip.label)],1)
    while (replace%in%dontchoose) replace<-sample(tree$tip.label[grepl(genus,tree$tip.label)],1)
    tree$tip.label[tree$tip.label==replace]<-prob_species[i]
  }
  return(tree)
}

load_zanne_tree<-function(){
  
  require(phytools)
  require(diversitree)
  
  zae <- read.tree("data/zanne1.1.tre")
  stemSamples <- load_stemSamples()
  wood_names<-sub(" ", "_", unique(stemSamples$Binomial))
  set.seed(42)
  zae_mod<-fix_problem_species(zae, 
                               prob_species = c("Melaleuca_decora","Petrophile_pulchella","Persoonia_nutans","Callistemon_linearis","Eucalyptus_sclerophylla","Isopogon_anemonifolius","Hakea_sericea","Olax_stricta","Jacksonia_scoparia"), 
                               dontchoose = wood_names)
  zanneTree<-diversitree:::drop.tip.fixed(phy = zae_mod,zae_mod$tip.label[!zae_mod$tip.label%in%wood_names])
  
  return(zanneTree)
  
}



#---------------------------------------------------------#
# load wood trait data

load_waterPercent.perGwetmass<-function(){
  
  require(dplyr)
  
  # read in initial covariate data
  covar.big <-read.csv('data/covariates_bigStems.csv', stringsAsFactor = F)
  covar.small <-read.csv('data/covariates_smallStems.csv', stringsAsFactor = F)
  
  # (1) calculate water content for each species, size class
  # g water per g dry mass
  water.percent <- c(with(covar.big, ((Fresh.mass..g. - Dry.mass..g.) / Fresh.mass..g.) * 100),
                     with(covar.small, ((Fresh.mass..g. - Dry.mass.total..g.) / Fresh.mass..g.)*100))
  water.percent <- data.frame(code=c(covar.big$Species, covar.small$Species),
                              StemSize=factor(c(rep('large', nrow(covar.big)), rep('small', nrow(covar.small)))),
                              water.percent, stringsAsFactors=F)
  
  ## aggregate by code
  group_by(water.percent, code) %>%
    summarize(meanWaterPerc = mean(water.percent, na.rm=TRUE),
              sdWaterPerc = sd(water.percent, na.rm=TRUE)) -> waterPercent
  
  #remove sd cols
  waterPercent<-waterPercent[,c("code","meanWaterPerc")]
  colnames(waterPercent)<-c("code","waterperc")
  length(unique(waterPercent$code)) #34
  
  return(waterPercent)
  
}

load_densityNbarkthick<-function(){
  
  require(dplyr)
  
  # read in initial covariate data
  covar.small <-read.csv('data/covariates_smallStems.csv', stringsAsFactor = F)
  
  # calculate wood density and bark thickness for each species
  #only have these values measured on small stems
  
  covar.small$density.gpercm3 <- with(covar.small, Dry.mass.wood..g. / Volume..g.)
  covar.small$barkthickness.mm <- with(covar.small, Diameter.wbark..mm. - Diameter.nobark..mm.)
  
  ## aggregate by species (all small size)
  group_by(covar.small, Species) %>%
    summarize(meanDensity = mean(density.gpercm3, na.rm=TRUE),
              sdDensity = sd(density.gpercm3, na.rm=TRUE),
              meanBarkthick = mean(barkthickness.mm, na.rm=TRUE),
              sdBarkthick = sd(barkthickness.mm, na.rm=TRUE)) -> densityNbarkthick
  colnames(densityNbarkthick)[1]<-"code"
  
  #remove sd cols
  densityNbarkthick<-densityNbarkthick[,c("code","meanDensity","meanBarkthick")]
  colnames(densityNbarkthick)<-c("code","density","barkthick")
  length(unique(densityNbarkthick$code)) #22
  
  return(densityNbarkthick)
  
}

load_XRF<-function(){
  
  require(dplyr)
  
  data <- read.csv('data/sequencing_T0/NextGenSeqencing_2016_Sample_MasterSpreadsheet.csv', stringsAsFactors=F)
  data <- data[!data$SampleCode == 'blank', ] # extract 'blank' and 'mock' samples from 'mat.otu', delete from 'data'
  
  # create dataframe containing metadata for initial sequencing and XRF data
  meta <- data[, c('SampleCode', 'StemSize', 'mgSample', 'NucleicAcidConc', 'ExtractionDate')]
  meta$code <- substr(meta$SampleCode, 1, 4)
  # add column to 'meta' indicating whether data obtained from independent/composite sample
  meta$compositeSample <- T
  meta$compositeSample[grep('[0-9A-Z]$', meta$SampleCode)] <- F #if there is a number on the back of the code, then it was composited
  
  #isolate XRF cols
  df.xrf<-data.frame(SampleCode=data$SampleCode, data[, 26:ncol(data)])
  indx<-meta[,c("SampleCode","StemSize","code","compositeSample")]
  df.xrf1<-left_join(indx,df.xrf)
  
  group_by(df.xrf1, code) %>%
    summarize(compos=paste(unique(compositeSample), collapse="_")) -> summ
  #why are there some species+size that are composited == TRUE and FALSE?
  
  ## aggregate by code
  group_by(df.xrf1, code) %>%
    summarize(meanP = mean(P, na.rm=TRUE),
              sdP = sd(P, na.rm=TRUE),
              
              meanK = mean(K, na.rm=TRUE),
              sdK = sd(K, na.rm=TRUE),
              
              meanCa = mean(Ca, na.rm=TRUE),
              sdCa = sd(Ca, na.rm=TRUE),
              
              meanMn = mean(Mn, na.rm=TRUE),
              sdMn = sd(Mn, na.rm=TRUE),
              
              meanFe = mean(Fe, na.rm=TRUE),
              sdFe = sd(Fe, na.rm=TRUE),
              
              meanZn = mean(Zn, na.rm=TRUE),
              sdZn = sd(Zn, na.rm=TRUE)
              
    ) -> xrf
  
  xrf<-xrf[,c("code","meanP","meanK","meanCa","meanMn","meanFe","meanZn")]
  colnames(xrf)<-c("code","P","K","Ca","Mn","Fe","Zn")
  
  return(xrf)
}

load_CN<-function(){
  
  require(dplyr)
  
  # read in CN data
  cndata <- read.csv('data/CN/JEFF_POWELL_CN_DATA_DIVERSITY_ROT_AUG_2013.csv', stringsAsFactors=F)
  colnames(cndata)<-c("sampleID","comments","mass","n.perc","c.perc")
  
  #create meta from XRF meta data
  data <- read.csv('data/sequencing_T0/NextGenSeqencing_2016_Sample_MasterSpreadsheet.csv', stringsAsFactors=F)
  data <- data[!data$SampleCode == 'blank', ] # extract 'blank' and 'mock' samples from 'mat.otu', delete from 'data'
  meta <- data[, c('SampleCode', 'StemSize', 'StemCode', 'mgSample', 'NucleicAcidConc', 'ExtractionDate')]
  meta$code <- substr(meta$SampleCode, 1, 4)
  meta$Stem<-as.numeric(substr(meta$SampleCode, 5, 6)) # these numbers correspond to Stem, right? NAs are because character missing or character was a capital letter (not number)
  
  # add column to 'meta' indicating whether data obtained from independent/composite sample
  meta$compositeSample <- T
  meta$compositeSample[grep('[0-9A-Z]$', meta$SampleCode)] <- F #if there is a number on the back of the code, then it was composited
  
  ## fix code for composited samples
  indx<-meta[,c("SampleCode","StemSize","Stem","code","compositeSample")]
  temp<-merge(cndata, indx, by.x="sampleID",by.y="SampleCode", all.x=TRUE)
  temp[is.na(temp$compositeSample),"compositeSample"]<-TRUE
  x<-temp[temp$compositeSample==TRUE,"sampleID"]
  xx<-unlist(strsplit(x,"pooled"))
  xx<-tolower(xx)
  xx[xx=="cali "]<-"cali"
  xx[xx=="acelextra"]<-"acel"
  temp[temp$compositeSample==TRUE,"code"]<-substr(xx, 1,4)
  temp[temp$compositeSample==TRUE,"StemSize"]<-"small"
  
  ## make codeStem col
  temp$codeStem<-paste(temp$code, temp$Stem, sep="")
  
  # aggregate values by code
  group_by(temp, code) %>%
    summarize(meanN = mean(n.perc, na.rm=TRUE),
              sdN = sd(n.perc, na.rm=TRUE),
              
              meanC = mean(c.perc, na.rm=TRUE),
              sdC = sd(c.perc, na.rm=TRUE)
              
    ) -> temp.agg
  aggSamps<-temp.agg[,c("code","meanN","meanC")]
  colnames(aggSamps)<-c("code","n.perc","c.perc")
  cn<-data.frame(aggSamps)
  length(unique(cn$code)) #33
  
  return(cn)
  
}

mergeTraitData<-function(){
  
  #load code-aggregated trait data
  waterPercent<-load_waterPercent.perGwetmass() ##### this is in units of g water per g of wet mass x 100
  densityNbarkthick<-load_densityNbarkthick()
  xrf<-load_XRF()
  cn<-load_CN()
  
  # merge together water percent, traits, and xrf.mean by 'code'
  tmp<-left_join(waterPercent, densityNbarkthick)
  tmp<-left_join(tmp, xrf)
  traits<-left_join(tmp, cn)
  
  # rename columns
  traits<-rename(traits, "N"="n.perc", "C"="c.perc")
  
  # use species-level small-stem estimates of density and barkthick for large-stem samples
  traits$species<-tolower(traits$code)
  traits$size<-"small"
  traits[tolower(traits$code)!=traits$code,"size"]<-"large"
  largeCodes<-as.data.frame(traits[traits$size=="large","code"])[,1]
  for(i in 1:length(largeCodes)){
    curr.code<-largeCodes[i]
    curr.species<-tolower(curr.code)
    filter(traits, species==curr.species) %>%
      filter(size=="small") -> curr.row
    curr.data<-data.frame(curr.row[,c("density","barkthick")])
    traits[traits$code == curr.code, c("density","barkthick")]<-curr.data
  }
  
  traits.code<-traits[,c("code", "species","size",
                         "waterperc","density","barkthick","P" ,"K" ,"Ca" ,"Mn","Fe","Zn" ,"N","C")]
  
  # #summarize trait ranges
  # traits.long<-as.data.frame(gather(traits, key=trait, value=value, -(1:3)))
  # group_by(traits.long, trait) %>%
  #   summarize(max=range(value, na.rm=TRUE)[1],
  #             min=range(value, na.rm=TRUE)[2])
  
  return(traits.code)
  
}



#---------------------------------------------------------#
# make a complete data subset

trim_OTUs <- function(minPerc, mat.otu){
  
  # trim OTUs that that rarely show up
  # define the minimum % of samples that an OTU needs to show up in
  # minPerc
  
  # calculate the minimum number of samples that an OTU needs to show up in
  numSamps <- dim(mat.otu)[1]
  minSamps <- floor( numSamps * (minPerc/100) )
  
  # remove OTUs that don't meet this criteria
  x <- apply(mat.otu > 0, 2, sum)
  mat.otu.trimmed <- mat.otu[, x > minSamps]
  
  return(mat.otu.trimmed)
  
}

subset_and_make_list<-function(seqSamples, traits.code, mat.otu){
  
  minPerc <- 20
  
  #annotate seqSamples with site and traits
  seqSamples %>%
    left_join(traits.code) %>%
    arrange(seq_sampName) %>%
    select(seq_sampName, code, species, site, size, waterperc, density, barkthick, C, N, P, K, Ca, Fe, Mn, Zn) -> mat.covars
  
  #make sure that the row order matches with mat.covars and mat.otu
  o<-match(mat.covars$seq_sampName, row.names(mat.otu))
  mat.otu.o<-mat.otu[o,]
  
  #get rid of rows with missing values in the covariates
  cov.cols<-c("species","site","size",
              "waterperc","density","barkthick",
              "C","N","P","K","Ca","Fe","Mn","Zn")
  mat.covars %>%
    filter(complete.cases(mat.covars[,cov.cols])) -> mat.covars.complete
  mat.otu.o[complete.cases(mat.covars[,cov.cols]),] -> mat.otu.complete
  sum(row.names(mat.otu.complete) != mat.covars.complete$seq_sampName) #this still needs to be 0
  
  #trim OTUs that do not show up in X percent of the samples
  mat.otu.trimmed <- trim_OTUs(minPerc = minPerc, mat.otu = mat.otu.complete)
  
  #save everything
  complete_subset_list<-list(otus = mat.otu.complete,
                             otus.trimmed = mat.otu.trimmed,
                             covariates = mat.covars.complete)
  
  if (any(
    mat.covars.complete$seq_sampName != row.names(mat.otu.complete)
    )) message("STOP: covariate data not matching otu matrix!")
  
  return(complete_subset_list)
  
}

save_complete_subset_list <- function(complete_subset_list){
  
  saveRDS(complete_subset_list, file = "derived_data/complete_subset_list.RData")
  
}



