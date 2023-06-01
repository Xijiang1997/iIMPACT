Installation
==========================

iIMPACT provides a novel method for spatial domain identification in spatially resolved transcriptomics (SRT) data, which integrates image and molecular profiles to improve the domain identification accuracy. It also has the ability to detect domain-specific spatially variable genes via a negative binomial regression model.

To implement iIMPACT function, you need to make sure that the following packages are already installed. ::

    library(ggplot2)
    library(Rcpp)
    library(RcppArmadillo)
    library(RcppDist)
    library(DirichletReg)
    library(mvtnorm)
    library(LaplacesDemon)
    library(SingleCellExperiment)
    library(scater)
    library(scran)

Load functions used to run iIMPACT.
::
    source('R/iIMPACT_functions.R')

