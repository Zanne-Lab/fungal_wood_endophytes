# concept figure
#######

#---------------------------------------------------------#
# make relative abundance table

make_commMats<-function(){
  
  sampleNames<-c("wood sp. 1, sample 1",
                 "wood sp. 1, sample 2",
                 "wood sp. 2, sample 1",
                 "wood sp. 2, sample 2")
  otuNames<-c("OTU a","OTU b")
  otumat<-mat.or.vec(nr=length(sampleNames), nc=length(otuNames))
  row.names(otumat)<-sampleNames
  colnames(otumat)<-otuNames
  
  # environmental filtering
  otu.modul1<-otumat
  otu.modul1["wood sp. 1, sample 1",]<-c(1, 0)
  otu.modul1["wood sp. 1, sample 2",]<-c(1, 0)
  otu.modul1["wood sp. 2, sample 1",]<-c(0, 2)
  otu.modul1["wood sp. 2, sample 2",]<-c(0, 2)
  
  # true checkerboard
  otu.modul2<-otumat
  otu.modul2["wood sp. 1, sample 1",]<-c(1, 2)
  otu.modul2["wood sp. 1, sample 2",]<-c(2, 1)
  otu.modul2["wood sp. 2, sample 1",]<-c(0, 0)
  otu.modul2["wood sp. 2, sample 2",]<-c(0, 0)
  
  # true anti-checkerboard
  otu.modul3<-otumat
  otu.modul3["wood sp. 1, sample 1",]<-c(1, 1)
  otu.modul3["wood sp. 1, sample 2",]<-c(2, 2)
  otu.modul3["wood sp. 2, sample 1",]<-c(0, 0)
  otu.modul3["wood sp. 2, sample 2",]<-c(0, 0)
  
  
  otu.list<-list(otu.modul1, otu.modul2, otu.modul3)
  names(otu.list)<-c("envfilter","checker","antichecker")
  return(otu.list)
  
}

plot_relAbund<-function(otumat){
  
  otumat.df<-data.frame(sampName=row.names(otumat), otumat)
  otumat.df.l<-gather(otumat.df, OTU, value, -sampName)
  otumat.df.l$sampName<-factor(otumat.df.l$sampName, levels=rev(levels(otumat.df.l$sampName)))
  otumat.df.l$OTU<-gsub(".", " ", otumat.df.l$OTU, fixed=TRUE)
  
  p<-ggplot(otumat.df.l, aes(x=OTU,y=sampName)) + 
    geom_tile(aes(fill = value), colour = "black") +
    geom_text(aes(label = value), color="black") +
    scale_fill_distiller(name="Rel. abundance", direction=1, palette="Purples") +
    labs(y = "", x = "") + 
    scale_x_discrete(expand = c(0, 0), position = 'top') +
    scale_y_discrete(expand = c(0, 0)) + 
    theme_linedraw() + theme(axis.ticks=element_blank(), plot.title=element_text(hjust=0.5, face="italic")) +
    guides(fill=FALSE)

  return(p)
}

#---------------------------------------------------------#
# make correlation tables

make_corsubs<-function(cor.list){
  
  # environmental filter
  corsub.modul1 <- cor.list[[1]]
  corsub.modul1['OTU a','OTU b'] <- -1 # shared env
  corsub.modul1['OTU b','OTU a'] <- 0 # res
  
  # true checker
  corsub.modul2 <- cor.list[[2]]
  corsub.modul2['OTU a','OTU b'] <- 1 # shared env
  corsub.modul2['OTU b','OTU a'] <- -1 # res
  
  # true anti-checker
  corsub.modul3 <- cor.list[[3]]
  corsub.modul3['OTU a','OTU b'] <- 1 # shared env
  corsub.modul3['OTU b','OTU a'] <- 1 # res
  
  corsub.list<-list(corsub.modul1, corsub.modul2, corsub.modul3)
  names(corsub.list)<-names(cor.list)
  return(corsub.list)
  
}

# extract_uniquePairDists in load_data.R (upper.tri)

extract_uniquePairDists_lowertri<-function(dist.mat){
  
  x<-dist.mat
  rowCol <- expand.grid(rownames(x), colnames(x))
  labs <- rowCol[as.vector(lower.tri(x,diag=F)),]
  df <- cbind(labs, x[lower.tri(x,diag=F)])
  colnames(df) <- c("sp1","sp2","dist")
  
  return(df)
}

plot_cor<-function(cormat){
  
  # take just the upper tri and make it long
  mat.df<-extract_uniquePairDists(cormat)
  colnames(mat.df)<-c("otu1","otu2","value")
  otuNames <- colnames(cormat)
  oneTOone<-data.frame(otu1 = otuNames, otu2 = otuNames, value = NA) #add back the 1:1 values
  mat.df<-rbind(mat.df, oneTOone)
  
  #fix plotting things
  mat.df$otu2<-factor(mat.df$otu2, levels=rev(otuNames))
  mat.df$value<-round(mat.df$value, digits=4)
  mat.df$sign <- " "
  mat.df[mat.df$value > 0 & !is.na(mat.df$value),"sign"]<-"+"
  mat.df[mat.df$value < 0 & !is.na(mat.df$value),"sign"]<-"-"
  mat.df$sign <- factor(mat.df$sign, levels = c("+","-"," "))
  
  #colors
  colorVec <- cor_colors()
  curr.colors <- c("-" = colorVec[['negRed']], 
                   "+" = colorVec[['posBlue']],
                   " " = "#bcbcbc")
  naGray<-"#bcbcbc"
  
  p<-ggplot(mat.df, aes(x=otu1,y=otu2)) + 
    geom_tile(aes(fill = sign), colour = "black") +
    geom_text(aes(label = sign), size=5) +
    labs(y = "", x = "") + 
    scale_x_discrete(expand = c(0, 0), position = 'top') +
    scale_y_discrete(expand = c(0, 0)) + 
    theme_linedraw() + 
    theme(axis.ticks=element_blank(), plot.title=element_text(hjust=0.5, face="italic"),
          panel.grid = element_blank()) +
    scale_fill_manual(values = curr.colors, na.value = naGray) +
    guides(fill=F) 
  p
  
  return(p)
}

plot_cor2<-function(cormat){
  
  cormat
  
  # take just the upper tri and make it long
  sharenv.df <- extract_uniquePairDists(cormat)
  sharenv.df
  colnames(sharenv.df) <- c("otu1","otu2","value")
  sharenv.df$corType <- "sharedEnv"
  
  # take just the lower tri and make it long
  rescor.df <- extract_uniquePairDists_lowertri(cormat)
  rescor.df
  colnames(rescor.df) <- c("otu1","otu2","value")
  rescor.df$corType <- "residual"
  
  #put everything together
  otuNames <- colnames(cormat)
  oneTOone<-data.frame(otu1=otuNames, otu2=otuNames, value=NA, corType="none") #add back the 1:1 values
  mat.df<-rbind(sharenv.df, rescor.df, oneTOone)
  
  #fix plotting things
  mat.df$otu2<-factor(mat.df$otu2, levels=rev(otuNames))
  mat.df$value<-round(mat.df$value, digits=4)
  mat.df$sign<- " "
  mat.df[mat.df$value > 0 & !is.na(mat.df$value),"sign"]<-"+"
  mat.df[mat.df$value < 0 & !is.na(mat.df$value),"sign"]<-"-"
  mat.df[mat.df$value == 0 & !is.na(mat.df$value),"sign"]<-"0"
  mat.df$sign <- factor(mat.df$sign, levels = c("+","-","0"," "))
  
  #colors
  colorVec <- cor_colors()
  curr.colors <- c("-" = colorVec[['negRed']], 
                   "+" = colorVec[['posBlue']],
                   " " = "#bcbcbc",
                   "0" = colorVec[['zeroGray']])
  naGray<-"#bcbcbc"
  
  p<-ggplot(mat.df, aes(x=otu1,y=otu2)) + 
    geom_tile(aes(fill = sign), colour = "black") +
    geom_text(aes(label=sign), size=5) +
    labs(y = "", x = "") + 
    scale_x_discrete(expand = c(0, 0), position = 'top') +
    scale_y_discrete(expand = c(0, 0)) + 
    theme_linedraw() + 
    theme(axis.ticks=element_blank(), plot.title=element_text(hjust=0.5, face="italic"),
          panel.grid = element_blank()) +
    scale_fill_manual(values = curr.colors, na.value = naGray) +
    guides(fill=FALSE) 
  
  return(p)
}


#---------------------------------------------------------#
# make correlation distribution figures

make_corrdf<-function(mean, sd, signifLim){
  
  lowerlim <- signifLim *-1
  upperlim <- signifLim
  
  numPairs <- 1000
  corr <- rnorm(numPairs, mean=mean, sd=sd)
  corr<-corr[corr <= 1 & corr >= -1] #get rid of 'correlation' values that are greater than 1 or less than -1
  corr.df<-data.frame(corr)
  
  corr.df$signif<-"zero"
  corr.df[corr.df$corr >= upperlim,"signif"]<-"pos"
  corr.df[corr.df$corr <= lowerlim,"signif"]<-"neg"
  
  return(corr.df)
  
}

make_distbdfs<-function(){
  
  signifLim<-0.2
  
  #environmental filter
  corr.df1<-make_corrdf(mean= -.5, sd=.2, signifLim=signifLim)
  env.df1<-make_corrdf(mean=-.5, sd=.2, signifLim=signifLim)
  res.df1<-make_corrdf(mean=0, sd=.1, signifLim=signifLim)
  
  #true checkerboard
  corr.df2<-make_corrdf(mean = .6, sd=.2, signifLim=signifLim)
  env.df2<-make_corrdf(mean = .5, sd=.2, signifLim=signifLim)
  res.df2<-make_corrdf(mean = -.3, sd=.2, signifLim=signifLim)
  
  #true anti-checkerboard
  corr.df3<-make_corrdf(mean = .6, sd=.2, signifLim=signifLim)
  env.df3<-make_corrdf(mean = .5, sd=.2, signifLim=signifLim)
  res.df3<-make_corrdf(mean = .3, sd=.2, signifLim=signifLim)
  
  corr.df.list<-list(corr.df1, corr.df2, corr.df3)
  env.df.list<-list(env.df1, env.df2, env.df3)
  res.df.list<-list(res.df1, res.df2, res.df3)
  
  names(corr.df.list)<-names(env.df.list)<-names(res.df.list)<-c("envfilter","checker","antichecker")
  distrib.dfs<-list(corr=corr.df.list, env=env.df.list, res=res.df.list)
  return(distrib.dfs)
}

concept_theme<-function(){
  
  require(ggplot2)
  
  concept.theme <- theme_bw(base_size = 10, base_family = "Helvetica") +
    theme(panel.border = element_rect(colour = "black"),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text.x = element_text(face='bold.italic', hjust=0.05),
          strip.text.y = element_text(face='bold.italic', hjust=0.05),
          strip.background = element_rect(fill = 'white', colour='black'),
          legend.key = element_blank(),
          plot.title=element_text(hjust=0, vjust=0.5, face='bold'),
          plot.margin = unit(c(0.05,0.05,0.05,0.05),"in"),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.x = element_text(size = 18)
    )
  
  return(concept.theme)
  
}

plot_cor_distrib<-function(corr.df, xlabel){
  
  #colors
  colorVec <- cor_colors()
  curr.colors <- c("neg" = colorVec[['negRed']], 
                   "pos" = colorVec[['posBlue']],
                   "zero" = colorVec[['zeroGray']])
  
  #plots
  mytheme<-concept_theme()
  p.corr<-ggplot(corr.df, aes(x=corr, fill=signif)) + 
    geom_histogram() +
    geom_vline(xintercept = 0, linetype=2) +
    xlab(xlabel) + ylab(" ") +
    scale_fill_manual(values = curr.colors) + 
    mytheme + guides(fill=FALSE) +
    scale_x_continuous(breaks = c(-1,0,1), limits = c(-1,1))
  
  return(p.corr)
  
}

