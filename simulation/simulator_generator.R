generate_natural_mat <- function(cell_pop, gene_pop, n_each, d_each, sigma, modifier, tol = 1){
  h <- nrow(cell_pop)
  cell_mat <- do.call(rbind, lapply(1:h, function(x){
    pos <- stats::runif(n_each)
    cbind(pos*cell_pop[x,1] + (1-pos)*cell_pop[x,3] + stats::rnorm(n_each, sd = sigma),
          pos*cell_pop[x,2] + (1-pos)*cell_pop[x,4] + stats::rnorm(n_each, sd = sigma))
  }))
  cell_mat <- pmax(cell_mat, tol)
  n <- nrow(cell_mat)
  k <- ncol(cell_mat)

  # construct the gene information
  g <- nrow(gene_pop)
  gene_mat <- do.call(rbind, lapply(1:g, function(x){
    pos <- stats::runif(d_each)
    cbind(pos*gene_pop[x,1] + (1-pos)*gene_pop[x,3] + stats::rnorm(d_each, sd = sigma),
          pos*gene_pop[x,2] + (1-pos)*gene_pop[x,4] + stats::rnorm(d_each, sd = sigma))
  }))
  gene_mat <- pmax(gene_mat, tol)
  d <- nrow(gene_mat)

  res <- .reparameterize(cell_mat*sqrt(modifier), gene_mat*sqrt(modifier))

  list(nat_mat = res$u_mat %*% t(res$v_mat), cell_mat =  res$u_mat, gene_mat = res$v_mat)
}

#####################

# observe that for generator_pcmf_poisson, nat_mat isn't really the natural parameters -- it's instead the matrix of mean parameters themselves
generator_pcmf_poisson <- function(nat_mat, dropout_prob = 0.1, ...){
  stopifnot(all(nat_mat > 0))
  
  n <- nrow(nat_mat); d <- ncol(nat_mat)

  obs_mat <- matrix(0, ncol = d, nrow = n)
  dropout_mat <- matrix(0, ncol = d, nrow = n)

  for(i in 1:n){
    for(j in 1:d){
      obs_mat[i,j] <- stats::rpois(1, nat_mat[i,j])
      dropout_mat[i,j] <- stats::rbinom(1, size = 1, prob = 1-dropout_prob)
    }
  }

  # 1 means not dropped
  obs_mat[dropout_mat == 0] <- 0

  obs_mat
}


generator_esvd_poisson <- function(nat_mat,...){
  n <- nrow(nat_mat); d <- ncol(nat_mat)

  obs_mat <- matrix(0, ncol = d, nrow = n)
  for(i in 1:n){
    for(j in 1:d){
      obs_mat[i,j] <- stats::rpois(1, exp(nat_mat[i,j]))
    }
  }

  obs_mat
}

##########

# draw a negative binomial from the exponential of inner product
generator_zinb_nb <- function(nat_mat, r_vec = rep(100, ncol(nat_mat)), ...){
  stopifnot(all(nat_mat < 0))
  
  n <- nrow(nat_mat); d <- ncol(nat_mat)

  obs_mat <- matrix(0, ncol = d, nrow = n)
  dropout_mat <- matrix(0, ncol = d, nrow = n)
  for(i in 1:n){
    for(j in 1:d){
      p <- exp(nat_mat[i,j])/(exp(nat_mat[i,j])+r_vec[j])
      r <- r_vec[j]
      obs_mat[i,j] <- stats::rnbinom(1, size = r, prob = p)

      dropout_mat[i,j] <- stats::rbinom(1, size = 1, prob = 1/(1+exp(-nat_mat[i,j])))
    }
  }

  # 1 means not dropped
  obs_mat[dropout_mat == 0] <- 0

  obs_mat
}

generator_esvd_nb <- function(nat_mat, scalar = 100,  ...){
  stopifnot(all(nat_mat < 0))
  
  n <- nrow(nat_mat); d <- ncol(nat_mat)

  obs_mat <- matrix(0, ncol = d, nrow = n)
  for(i in 1:n){
    for(j in 1:d){
      p <- 1-exp(-nat_mat[i,j]) # remember the p param in rnbinom is "flipped" in R
      obs_mat[i,j] <- stats::rnbinom(1, size = scalar, prob = p)
    }
  }

  obs_mat
}

#############

# draw an exponential from the negative of the natural parameter
generator_exponential <- function(nat_mat, ...){
  stopifnot(all(nat_mat > 0))
  
  n <- nrow(nat_mat); d <- ncol(nat_mat)

  obs_mat <- matrix(0, ncol = d, nrow = n)
  for(i in 1:n){
    for(j in 1:d){
      obs_mat[i,j] <- stats::rexp(1, rate = nat_mat[i,j])
    }
  }

  obs_mat
}

# draw a gaussian from the natural parameter
generator_gaussian <- function(nat_mat, sd_val = 0.25, tol = 1e-3, ...){
  n <- nrow(nat_mat); d <- ncol(nat_mat)

  obs_mat <- matrix(0, ncol = d, nrow = n)
  for(i in 1:n){
    for(j in 1:d){
      obs_mat[i,j] <- stats::rnorm(1, nat_mat[i,j], sd = sd_val)
    }
  }

  obs_mat[obs_mat < 0] <- tol

  obs_mat
}

# draw a gaussian from the negative of the natural parameter
generator_curved_gaussian <- function(nat_mat, scalar = 2, tol = 1e-3, ...){
  stopifnot(all(nat_mat > 0))
  n <- nrow(nat_mat); d <- ncol(nat_mat)

  obs_mat <- matrix(0, ncol = d, nrow = n)
  for(i in 1:n){
    for(j in 1:d){
      obs_mat[i,j] <- stats::rnorm(1, 1/nat_mat[i,j], sd = 1/(scalar*nat_mat[i,j]))
    }
  }

  obs_mat[obs_mat < 0] <- tol

  obs_mat
}

#################################################

.identification <- function(cov_x, cov_y, check = F, tol = 1e-6){
  stopifnot(all(dim(cov_x) == dim(cov_y)), nrow(cov_x) == ncol(cov_x))
  if(nrow(cov_x) == 1){
    return(matrix((as.numeric(cov_y)/as.numeric(cov_x))^(1/4), 1, 1))
  }
  
  eigen_x <- eigen(cov_x)
  eigen_y <- eigen(cov_y)
  
  Vx <- eigen_x$vectors
  Vy <- eigen_y$vectors
  
  if(any(eigen_x$values <= tol) | any(eigen_y$values <= tol)) warning("Detecting rank defficiency in reparameterization step")
  
  Dx <- diag(eigen_x$values)
  Dy <- diag(eigen_y$values)
  
  tmp <- sqrt(Dy) %*% t(Vy) %*% Vx %*% sqrt(Dx)
  svd_tmp <- svd(tmp)
  R <- svd_tmp$u %*% t(svd_tmp$v)
  
  # run a check
  if(check){
    Q <- t(R) %*% tmp
    stopifnot(sum(abs(Q - t(Q))) <= 1e-6)
  }
  
  Dx_inv <- Dx; diag(Dx_inv) <- 1/diag(Dx)
  sym_prod <- Vx %*% sqrt(Dx_inv) %*% t(R) %*% sqrt(Dy) %*% t(Vy)
  sym_prod[which(abs(sym_prod) <= tol)] <- 0
  
  if(check){
    stopifnot(sum(abs(sym_prod - t(sym_prod))) <= 1e-6)
  }
  
  eigen_sym <- eigen(sym_prod)
  T_mat <- diag(sqrt(eigen_sym$values)) %*% t(eigen_sym$vectors)
  
  if(check){
    mat1 <- T_mat %*% cov_x %*% t(T_mat)
    
    T_mat_inv <- solve(T_mat)
    mat2 <- t(T_mat_inv) %*% cov_y %*% T_mat_inv
    stopifnot(sum(abs(mat1 - mat2)) <= 1e-6)
  }
  
  # adjust the transformation so it yields a diagonal matrix
  eig_res <- eigen(T_mat %*% cov_x %*% t(T_mat))
  
  t(eig_res$vectors) %*% T_mat
}


#' Function to reparameterize two matrices
#'
#' test
#'
#' Designed to output matrices of the same dimension as \code{u_mat}
#' and \code{v_mat}, but linearly transformed so \code{u_mat \%*\% t(v_mat)}
#' is preserved but either \code{u_mat \%*\% t(u_mat)} is diagonal and equal to
#' \code{v_mat \%*\% t(v_mat)} (if \code{equal_covariance} is \code{FALSE})
#' or \code{u_mat \%*\% t(u_mat)/nrow(u_mat)} is diagonal and equal to
#' \code{v_mat \%*\% t(v_mat)/nrow(v_mat)} (if \code{equal_covariance} is \code{TRUE})
#'
#' @param u_mat matrix of dimension \code{n} by \code{k}
#' @param v_mat matrix of dimension \code{p} by \code{k}
#' @param equal_covariance boolean
#'
#' @return list of two matrices
.reparameterize <- function(u_mat, v_mat, equal_covariance = F){
  stopifnot(ncol(u_mat) == ncol(v_mat))
  n <- nrow(u_mat); p <- nrow(v_mat)
  
  res <- .identification(t(u_mat) %*% u_mat, t(v_mat) %*% v_mat)
  
  if(equal_covariance){
    list(u_mat = (n/p)^(1/4)*u_mat %*% t(res), v_mat = (p/n)^(1/4)*v_mat %*% solve(res))
  } else {
    list(u_mat = u_mat %*% t(res), v_mat = v_mat %*% solve(res))
  }
  
}



