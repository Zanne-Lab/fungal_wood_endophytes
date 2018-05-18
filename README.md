Fungal endophyte communities rarely exhibit competitive exclusion patterns across a range of wood habitats
================

To determine how endophytic fungal taxa respond to both abiotic and biotic environmental drivers, we characterized these communities using high-throughput amplicon sequencing and measured wood functional traits in stems sampled from 22 species of plants growing in woodlands near Richmond, NSW, Australia.

This repository holds the data and code needed to reproduce our analyses, but users should be aware that two types of analyses, (i) latent variable co-occurrance models run with the R package boral and (ii) model-based classification run with the R package RCPmod are quite computationally-intensive (e.g. >24 hrs processing time). Model fitting for these analyses were conducted with the help of cluster computing. Scripts used to do so are included in the folder `forCluster`.  Intermediate R objects from those processes are included in the `derived_data` folder and are incorporated into the analysis workflows described below.

### Re-running the workflow...

Users can re-run our analyses by downloading this repository to their local computer and either running (1) a process mediated by the R package `remake` or (2) an old-fashioned R script with the same content. The `remake` package is a workflow manager can be very useful for analyses that are complex because it does not re-run processes for which there is already a cached data product or intermediate.  

The main inputs are the data associated with the harvest of the wood, both the meta-genomics data and the characterization of the wood substrates.  Another input for an analysis is the plant phylogeny previously published as Zanne et al. (2015).  Each step in the process is an R function that can be found in the `code` folder. The intermediate data products, e.g from cluster analyses, are located in the `derived_data` folder. Output objects include tables and figures that can be found in the folder "output".  

We recommend working in `RStudio`https://www.rstudio.com/products/rstudio/download/. If using `RStudio`, opening the `.Rproj` file will automatically set the working directory to the top folder, ensuring that paths to files will work.


#### Option 1. Use the remake package, which is "Make-like build management"

To re-run this workflow first download the R package called `remake` from https://github.com/richfitz/remake using this command:
``` r
install.packages("devtools")
devtools::install_github("richfitz/remake")
```

Next, make sure that you have all of the required R packages that are listed at the top of the yaml file `remake.yml`. Here is a one-liner for that
```r
remake::install_missing_packages()
```

Now, you should be able to re-run all of the analyses in `R` using this command:
```r
remake::make()
```
For a remake tutorial, go to https://github.com/ropenscilabs/remake-tutorial


#### Option 2. Run the "linear_remake.R" script

This is a traditional R script that uses the same R functions and data files as in option 1, above.  Please be aware that some of the processes in this script can take a while to run.





