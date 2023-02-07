Installation
==========================

iIMPACT provides a novel method for spatial domain identification in spatial transcriptomics (ST) data, which integrates image and molecular profiles to improve the domain identification accuracy. 

To implement iIMPACT function, you must make sure that the following packages are already installed. ::

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

Load functions used to run iIMPACT.::

    source('R/iIMPACT_functions.R')

