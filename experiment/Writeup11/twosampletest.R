compute_twosample_pvalue <- function(factor_vec,
                                     mean_vec, var_vec,
                                     k = 3,
                                     verbose = 0){
  stopifnot(length(factor_vec) == length(mean_vec),
            length(factor_vec) == length(var_vec))
  n <- length(mean_vec)

  if(verbose > 0) print("Computing pairwise distances")
  dist_mat <- .compute_pairwise_distances(mean_vec = mean_vec,
                                          var_vec = var_vec)

  if(verbose > 0) print("Computing MST")
  edge_mat <- .get_mst(dist_mat = dist_mat, k = k)
  sparse_mat <- .convert_edge2dgCMatrix(edge_mat = edge_mat, n = n)

  if(verbose > 0) print("Computing graph statistics")
  stat_list <- .gather_graph_statistics(factor_vec = factor_vec,
                                        sparse_mat = sparse_mat)
  pval <- stats::pnorm(-stat_list[["Zw"]])

  list(test_stat = stat_list[["Zw"]], pval = pval)
}

## https://djalil.chafai.net/blog/2010/04/30/wasserstein-distance-between-two-gaussians/
.compute_pairwise_distances <- function(mean_vec, var_vec){
  stopifnot(length(mean_vec) == length(var_vec), all(var_vec > 0))

  n <- length(mean_vec)
  dist_mat <- matrix(0, n, n)
  for(i in 1:(n-1)){
    for(j in 2:n){
      term1 <- (mean_vec[i] - mean_vec[j])^2
      term2 <- (var_vec[i] + var_vec[j] - 2*(var_vec[i]*var_vec[j])^(1/2))
      dist_mat[i,j] <- sqrt(term1 + term2)
    }
  }

  dist_mat <- dist_mat + t(dist_mat)
  stats::as.dist(dist_mat)
}

# returns a 2-column matrix
.get_mst <- function(dist_mat, k){
  res <- ade4::mstree(xdist = dist_mat, ngmax = k)
  class(res) <- "matrix"
  attr(res,"degrees") <- NULL
  attr(res, "call") <- NULL
  res
}

.convert_edge2dgCMatrix <- function(edge_mat, n){
  i_vec <- edge_mat[,1]
  j_vec <- edge_mat[,2]
  sparse_mat <- Matrix::sparseMatrix(i = i_vec, j = j_vec,
                                     x = rep(1, length(i_vec)),
                                     dims = c(n,n))
  sparse_mat <- sparse_mat + Matrix::t(sparse_mat)
  sparse_mat@x <- rep(1, length(sparse_mat@x))
  sparse_mat
}

# using calculations in gTests package, https://cran.r-project.org/web/packages/gTests/gTests.pdf
# https://arxiv.org/pdf/1604.06515.pdf
.gather_graph_statistics <- function(factor_vec, sparse_mat){
  stopifnot(length(factor_vec) == nrow(sparse_mat),
            inherits(sparse_mat, "dgCMatrix"),
            is.factor(factor_vec), length(levels(factor_vec)) == 2)
  N <- nrow(sparse_mat)

  idx1 <- which(factor_vec == levels(factor_vec)[1])
  idx2 <- setdiff(1:N, idx1)

  n1 <- length(idx1)
  n2 <- length(idx2)

  R1 <- sum(sparse_mat[idx1, idx1])/2
  R2 <- sum(sparse_mat[idx2, idx2])/2

  nE <- sum(sparse_mat)/2
  node_deg_list <- lapply(1:N, function(i){
    .nonzero_col(sparse_mat, col_idx = i, bool_value = F)
  })
  total_deg <- sapply(node_deg_list, length)
  nEi <- sum(total_deg * (total_deg - 1)) # pair of nodes sharing a node * 2

  mu1 <- nE * (n1*(n1-1)) / (N*(N-1))
  mu2 <- nE * (n2*(n2-1)) / (N*(N-1))
  V1 <- nEi * (n1*(n1-1)*(n1-2)) / (N*(N-1)*(N-2)) + (nE*(nE-1)-nEi) * (n1*(n1-1)*(n1-2)*(n1-3)) / (N*(N-1)*(N-2)*(N-3)) + mu1 - mu1^2
  V2 <- nEi * (n2*(n2-1)*(n2-2)) / (N*(N-1)*(N-2)) + (nE*(nE-1)-nEi) * (n2*(n2-1)*(n2-2)*(n2-3)) / (N*(N-1)*(N-2)*(N-3)) + mu2 - mu2^2
  V12 <- (nE*(nE-1)-nEi) * (n1*n2*(n1-1)*(n2-1)) / (N*(N-1)*(N-2)*(N-3)) - mu1*mu2

  Zw <- (n2*(R1-mu1)+n1*(R2-mu2)) / sqrt(n2^2*V1 + n1^2*V2 + 2*n1*n2*V12)

  tmp <- list(R1 = R1, R2 = R2,
              n1 = n1, n2 = n2, N = N,
              nE = nE, nEi = nEi,
              mu1 = mu1, mu2 = mu2,
              V1 = V1, V2 = V2, V12 = V12,
              Zw = Zw)
  tmp
}

########################

.nonzero_col <- function(mat, col_idx, bool_value){
  stopifnot(inherits(mat, "dgCMatrix"), col_idx %% 1 == 0,
            col_idx > 0, col_idx <= ncol(mat))

  val1 <- mat@p[col_idx]
  val2 <- mat@p[col_idx+1]

  if(val1 == val2) return(numeric(0))
  if(bool_value){
    # return the value
    mat@x[(val1+1):val2]
  } else {
    # return the row index
    mat@i[(val1+1):val2]+1
  }
}
