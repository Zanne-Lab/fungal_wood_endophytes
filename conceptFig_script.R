# make concept figure

#load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(boral) #NOTE---- boral uses JAGS so you need to first install from http://sourceforge.net/projects/mcmc-jags/files/JAGS/

#load functions
source("code/conceptFig.R")
source("code/load_data.R")
source("code/plottingTheme.R")

# -------------------------------------------------------------------#
# conceptFig.R : Figure that conceptually links boral to checkerboard patterns


# relative abundance table
otu.list<-make_commMats()
pdf(file="output/conceptFig/relabund.pdf", width=8, height=3)
grid.arrange(
  plot_relAbund(otu.list[[1]]),
  plot_relAbund(otu.list[[2]]), 
  plot_relAbund(otu.list[[3]]), 
  ncol=3)
dev.off()

# correlation matrix
cor.list<-lapply(otu.list, cor)
pdf(file="output/conceptFig/cor.pdf", width=5, height=1.5)
grid.arrange(
  plot_cor(cormat=cor.list[[1]]),
  plot_cor(cormat=cor.list[[2]]),
  plot_cor(cormat=cor.list[[3]]), 
  ncol=3)
dev.off()

# partitioned correlation matrix
corsub.list<-make_corsubs(cor.list)
pdf(file="output/conceptFig/cor_partition.pdf", width=5, height=1.5)
grid.arrange(
  plot_cor2(cormat=corsub.list[[1]]), 
  plot_cor2(cormat=corsub.list[[2]]), 
  plot_cor2(cormat=corsub.list[[3]]), 
  ncol=3)
dev.off()

# distribution of correlations
distbdfs.list <- make_distbdfs()
xlabel<-"Shared environment correlation"
p.env1<-plot_cor_distrib(corr.df = distbdfs.list[['env']][[1]], xlabel=xlabel)
p.env2<-plot_cor_distrib(corr.df = distbdfs.list[['env']][[2]], xlabel=xlabel)
p.env3<-plot_cor_distrib(corr.df = distbdfs.list[['env']][[3]], xlabel=xlabel)
xlabel<-"Residual correlation"
p.res1<-plot_cor_distrib(corr.df = distbdfs.list[['res']][[1]], xlabel=xlabel)
p.res2<-plot_cor_distrib(corr.df = distbdfs.list[['res']][[2]], xlabel=xlabel)
p.res3<-plot_cor_distrib(corr.df = distbdfs.list[['res']][[3]], xlabel=xlabel)

pdf(file="output/conceptFig/distribcorrBoral.pdf", width=18, height=3)
grid.arrange(p.env1, p.res1,
             p.env2, p.res2, 
             p.env3, p.res3, ncol=6)
dev.off()

