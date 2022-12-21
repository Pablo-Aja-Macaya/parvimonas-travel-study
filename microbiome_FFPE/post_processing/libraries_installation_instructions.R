# Install libraries from new_filtersCCR89DEF.R file and analisis_R_datos_clinicos.R

# ---- INSTALAR BIBLIOTECAS ----
cran_packages <- c("bookdown", "knitr", "devtools", "igraph", "qgraph", 
                   "RColorBrewer", "tidyverse", "network", "sna", "grid", "gridExtra",
                   "Rcpp", "RcppArmadillo", "pulsar", "huge", "networkD3", "plyr", "remotes")
git_source <- c("zdk123/SpiecEasi", "hallucigenia-sparsa/seqtime", "briatte/ggnet", 
                "hadley/lazyeval", "hadley/dplyr", "jbisanz/qiime2R", "YuLab-SMU/ggtree",
                "kassambara/ggpubr", "microbiome/microbiome", "microsud/microbiomeutilities") # fuente/nombre del paquete
bio_packages <- c("phyloseq", "microbiome", "genefilter", "ggtree")
git_packages <- c("SpiecEasi", "seqtime", "ggnet") # nombre del paquete


# Instalar paquetes CRAN
.inst <- cran_packages %in% installed.packages()
if(any(!.inst)) {
  install.packages(cran_packages[!.inst])
}
# Intalar paquetes BioConductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
.inst <- bio_packages %in% installed.packages()
if(any(!.inst)) {
  BiocManager::install(bio_packages[!.inst])
}
# Instalar paquetes GitHub
.inst <- git_source %in% installed.packages()
if(any(!.inst)) {
  devtools::install_github(git_source[!.inst])
}
# #Cargar los paquetes instalados en lugar de con library
sapply(c(cran_packages, bio_packages, git_packages), require, character.only = TRUE)



### INDEPENDENT

## dplyr
# https://www.r-project.org/nosvn/pandoc/dplyr.html

# the latest released version from CRAN with

install.packages("dplyr")
# another option if devtools gives problems: the latest development version from github with

if (packageVersion("devtools") < 1.6) {
  install.packages("devtools")
}
devtools::install_github("hadley/lazyeval")
devtools::install_github("hadley/dplyr")



## phyloseq 
# https://www.bioconductor.org/packages/release/bioc/html/phyloseq.html

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("phyloseq")


## qiime2R 
# https://rdrr.io/github/jbisanz/qiime2R/
install.packages("remotes")
remotes::install_github("jbisanz/qiime2R")

# Bioconductor version
library(BiocManager)
BiocManager::install("ggtree")

# github version
devtools::install_github("YuLab-SMU/ggtree")

# Install
if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/ggpubr")
