#devtools and dependancies
install.packages("xml2") #install.packages("xml2", dependencies=FALSE, INSTALL_opts = c('--no-lock'))
install.packages("curl")
install.packages("rversions")
install.packages("openssl")
install.packages("httr")
install.packages("gh")
install.packages("covr")
install.packages("roxygen2")
install.packages("usethis")
install.packages("devtools")
library(devtools) 

#FishPlot
install_github("chrisamiller/fishplot")

#ClonEvol
install_github("hdng/clonevol")
install.packages("gridBase")
install.packages("gridExtra")
install.packages("ggplot2")
install.packages("igraph")
install.packages("packcircles")
install_github("hdng/trees")

#SciClone
#install IRanges from bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("IRanges")
BiocManager::install("limma")
install_github("genome/bmm")
install_github("genome/sciClone")