
# ***README***
# The following script is used to run the iIMPACT method for spatial transcriptomics data
# proposed in the submitted manuscript titled "XXXX"
# ***END***


# load cpp functions
Rcpp::sourceCpp('R/iIMPACT_cpp_function.cpp')

# load other packages
library(DirichletReg)
library(mvtnorm)
library(LaplacesDemon)
library(SingleCellExperiment)
library(scater)
library(scran)

# Summarize neighborhood information: Extract the indices of all neighbors for each spot
# loc: a n*2 matrix for x and y coordinates of spots, where n is the number of spots
# n_neighbor: number of neighbors for lattice 

get.neighbor <- function(loc, n_neighbor){
  P <- matrix(0, nrow = nrow(loc), ncol = n_neighbor)
  loc <- as.matrix(loc)
  if (n_neighbor == 4){
    loc <- round(loc)
    aa <- sqrt(2)
  } else if (n_neighbor == 6){
    aa <- sqrt(3)
  } else {aa <- 1.2}
  
  dist_matrix <- vectorized_pdist(loc, loc)
  min_dist <- min(dist_matrix[dist_matrix > 0])
  
  dist_threshold <- min_dist*(aa - 1)*0.5 + min_dist
  #print(min_dist)
  #print(dist_threshold)
  
  for (i in 1:nrow(loc)){
    k <- 1
    for (j in 1:nrow(loc)){
      if (dist_matrix[i, j] > 0 & dist_matrix[i, j] < dist_threshold){
        P[i, k] <- j
        k <- k + 1
      }}}
  return(P)
}

# Summarize cell abundance information: Extract the number of cells for different cell types for each spot
# cell_loc: a m*2 matrix for x and y coordinates of cells, where m is the number of cells
# cell_type: a vector with length m to record the cell types for all cells
# spot_loc: a n*2 matrix for x and y coordinates of spots
# lattice: lattice type, 'hexagon' or 'square'

get_cell_abundance <- function(cell_loc, cell_type, spot_loc, lattice = 'hexagon'){
  
  n_cell <- dim(cell_loc)[1]
  n_spot <- dim(spot_loc)[1]
  
  
  dist_spot <- vectorized_pdist(as.matrix(spot_loc), as.matrix(spot_loc))
  
  min_dist <- min(dist_spot[dist_spot > 0])
  
  if (lattice == 'hexagon'){
    threshold <- min_dist*1.01*sqrt(3)/2
  }
  else if(lattice == 'square'){
    threshold <- min_dist*1.01 / sqrt(2)
  }
  else {stop('Parameter lattice: Please enter a correct value.')}
  
  cell_names <- unique(cell_type)
  V <- matrix(0, ncol = length(cell_names) , nrow = n_spot)
  colnames(V) <- cell_names
  
  # assign cells to spots
  count <- 0
  for (i in 1:n_cell){
    ne_class <- which(cell_names == cell_type[i])
    dist_matrix <- vectorized_pdist(as.matrix(cell_loc[i,]), as.matrix(spot_loc))
    
    list_spot <- which(dist_matrix <= threshold)
    for (index_spot in list_spot){
      V[index_spot, ne_class] <- V[index_spot, ne_class] + 1
    }
    
    if (floor(i*100/n_cell) == count)
    {
      print(paste0(count, '% has been done'))
      count <- count + 10
    }
  }
  return(V)
}


# Process the gene expression count matrix: Normalize and get the low-dimensional representation of gene expression profile via PCA
# count: a n*p matrix for gene expression counts, where n is the number of spots and p is the number of genes
# n_PC: number of principle components for the gene expression profile
process_gene_expression <- function(count, n_PC = 3, n_HVG = 2000){
  n_gene <- ncol(count)
  n_spot <- nrow(count)
  
  if (is.null(colnames(count)) == FALSE){
  rowData <- data.frame(gene = colnames(count))} else {rowData <- data.frame(gene = 1:n_gene)}
  
  colnames(count) <- rowData$gene
  
  sce <- SingleCellExperiment(assays = list(counts = as(t(count), 'dgCMatrix')), rowData = rowData)
  
  sce <- logNormCounts(sce)
  
  n_HVG <- min(n_HVG, n_gene)
  
  set.seed(100)
  dec <- scran::modelGeneVar(sce)
  top <- scran::getTopHVGs(dec, n = n_HVG)
  
  set.seed(101)
  sce <- scater::runPCA(sce, subset_row = top, ncomponents = n_PC)
  
  Y <- sce@int_colData@listData$reducedDims@listData$PCA
  
  return(Y)
}

# run iIMPACT
# V: spot-level cell abundance (n*q matrix, q is the number of cell types)
# Y: low-dimensional representation of gene expression profile via PCA
# G: obtained from get.neighbor function
# n_cluster: number of spatial domains
# w: scaling parameter for image profile
# label_switch_refer, an index of column for label switching, usually choose column index of tumor cell type
run_iIMPACT <- function(V, Y, G, n_cluster, w, label_switch_refer = 1){
  n_spot <- nrow(Y)
  n_cell_type <- ncol(V)
  n_PC <- ncol(Y)
  set.seed(123)
  
  e <- rep(1, n_cluster)
  f <- rep(1, n_spot)
  omega_initial <- rep(1/n_cell_type, n_cell_type)
  alpha <- rep(1, n_cell_type)
  tau <- 0.01
  eta <- rep(0, n_PC)
  sigma_initial <- diag(1, ncol = n_PC, nrow = n_PC)
  mu_initial <- rmvnorm(1, eta, sigma_initial/tau)
  Z_initial <- matrix(sample(1:n_cluster, n_spot, replace = T), ncol = 1)
  alpha_gamma <- 0.1
  beta_gamma <- 0.1
  
  # run iIMPACT
  result <- iIMPACT(V, Y, G, n_cluster, e, f, eta, tau,  alpha_gamma, beta_gamma, mu_initial, diag(sigma_initial), alpha, omega_initial, Z_initial,rep(w, n_spot))
  
  Z_p <- result[["Z"]]
  Z_p_old <- Z_p
  omega_p <- result[['omega']]
  mean_p <- result[['mu']]
  
  n_iter <- nrow(mean_p)
  K <- n_cluster
  
  omega_p_new <- array(0, dim = c(floor(n_iter/2), K, n_cell_type))
  mean_p_new <- array(0, dim = c(floor(n_iter/2), K, n_PC))
  
  # switch label
  for (i in (n_iter - floor(n_iter/2) + 1):n_iter){
    mean_now <- omega_p[i,seq(label_switch_refer, (K - 1)*n_cell_type + label_switch_refer, n_cell_type)]
    for (j in 1:K){
      Z_p[i, which(Z_p_old[i, ] == order(mean_now)[j])] <- j
      old_index = order(mean_now)[j]
      new_index = j
      omega_p_new[i - (n_iter - floor(n_iter/2)), new_index, ] <- omega_p[i, ((old_index - 1)*n_cell_type + 1):(old_index*n_cell_type)]
      mean_p_new[i - (n_iter - floor(n_iter/2)), new_index, ] <- mean_p[i, ((old_index - 1)*n_PC + 1):(old_index*n_PC)]
    }
  }
  
  print('100% has been done')
  
  return(list(Z = Z_p[(n_iter - floor(n_iter/2) + 1):n_iter, ], omega = omega_p_new, cell_type = colnames(V), n_cluster = K))
}

# Get spatial domain identification results
get_spatial_domain <- function(result){
  Z_p <- result[['Z']]
  spatial_domain <- apply(Z_p, 2, getmode)
  return(spatial_domain)
}

# Get domain-level cell proportion
get_domian_cell_prop <- function(result){
  omega_p <- result[['omega']]
  domain_cell_proportion <- apply(omega_p, c(2, 3), median)
  colnames(domain_cell_proportion) <- result[['cell_type']]
  
  return(domain_cell_proportion)
}

# Get interactive zone
# prop_cut: threshold to define interactive zone
get_interactive_zone <- function(result, prop_cut = 0.9){
  Z_p <- result[['Z']]
  K <- result[['n_cluster']]
  Z_prop <- matrix(0, ncol = K, nrow = ncol(Z_p))
  for (i in 1:ncol(Z_p)){
    for (k in 1:K){
      Z_prop[i, k] <- mean(Z_p[, i] == k)
    }
  }
  
  interactive_zone <- apply(Z_prop, 1, max) < prop_cut
  
  return(interactive_zone)
}



                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 

# Helper functions
######################

vectorized_pdist <- function(A,B){
  an = apply(A, 1, function(rvec) crossprod(rvec,rvec))
  bn = apply(B, 1, function(rvec) crossprod(rvec,rvec))
  m = nrow(A)
  n = nrow(B)
  tmp = matrix(rep(an, n), nrow=m)
  tmp = tmp +  matrix(rep(bn, m), nrow=m, byrow=TRUE)
  sqrt( tmp - 2 * tcrossprod(A,B) )
}

getmode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
