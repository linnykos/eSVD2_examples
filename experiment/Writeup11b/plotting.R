histogram_plot <- function(col_template_vec,
                           gene_list,
                           teststat_vec,
                           bool_separate = F,
                           cex_legend = 0.6,
                           hist_spacing = 0.1,
                           legend_position = "topright",
                           max_abs_val = 30,
                           xlab = "Z-score",
                           ylab = "Frequency",
                           xlim = NULL,
                           ...){
  stopifnot(length(col_template_vec) == length(gene_list))

  idx_list <- lapply(gene_list, function(gene_vec){
    which(names(teststat_vec) %in% gene_vec)
  })
  col_vec <- rep(rgb(0.5, 0.5, 0.5, 0.5), length(teststat_vec))
  for(i in 1:length(idx_list)){
    col_vec[idx_list[[i]]] <- col_template_vec[i]
  }
  shuf_idx <- unlist(idx_list)
  shuf_idx <- sample(shuf_idx)

  teststat_vec <- pmax(pmin(teststat_vec, abs(max_abs_val)), -abs(max_abs_val))
  max_val <- max(abs(teststat_vec))
  if(all(is.null(xlim))) xlim <- c(-max_val, max_val)

  break_vec <- seq(-max_val-0.05, max_val+0.05, by = hist_spacing)
  break_vec <- c(break_vec, max(break_vec)+hist_spacing)

  if(!bool_separate){
    hist(teststat_vec, breaks = break_vec,
         xlim = xlim,
         xlab = xlab, ylab = ylab,
         freq = T, ...)
    lines(rep(0,2), c(0, 1e5), lwd = 1, lty = 3)
    for(i in shuf_idx){
      # note: needed since rug can not take in multiple colors at a time
      rug(teststat_vec[i], col = col_vec[i], lwd = 2)
    }

    legend(legend_position, names(gene_list),
           fill = col_template_vec, cex = cex_legend)
  } else {
    for(i in 1:length(gene_list)){
      idx <- idx_list[[i]]
      hist(teststat_vec[idx], breaks = break_vec,
           xlim = xlim,
           main = names(gene_list)[i],
           xlab = xlab, ylab = ylab,
           freq = T, ...)
      lines(rep(0,2), c(0, 1e5), lwd = 1, lty = 3)
      rug(teststat_vec[idx], col = col_vec[idx], lwd = 2)
    }
  }

  invisible()
}

plot_scatterplot_mean <- function(mat,
                                  esvd_res,
                                  nuisance_vec,
                                  case_control_variable,
                                  alpha_max = 50,
                                  asp = T,
                                  bool_logscale = F,
                                  cex = 1,
                                  col_points = rgb(0.5, 0.5, 0.5, 0.1),
                                  gene_names = NULL,
                                  max_num = 1e5,
                                  mean_type = "predicted",
                                  nuisance_lower_quantile = 0.01,
                                  only_nonzero = T,
                                  xlim = NA,
                                  ylim = NA,
                                  verbose = T,
                                  ...){
  stopifnot(mean_type %in% c("predicted", "posterior"))

  if(mean_type == "predicted"){
    nat_mat1 <- tcrossprod(esvd_res$x_mat, esvd_res$y_mat)
    nat_mat2 <- tcrossprod(esvd_res$covariates,
                           esvd_res$b_mat)
    mean_mat <- exp(nat_mat1 + nat_mat2)
  } else {
    res <- eSVD2:::compute_posterior(mat = mat,
                                     esvd_res = esvd_res,
                                     nuisance_vec = nuisance_vec,
                                     case_control_variable = case_control_variable,
                                     alpha_max = alpha_max,
                                     nuisance_lower_quantile = nuisance_lower_quantile)
    mean_mat <- res$posterior_mean_mat

    offset_var <- setdiff(colnames(esvd_res$covariates), case_control_variable)
    library_mat <- exp(tcrossprod(
      esvd_res$covariates[,offset_var],
      esvd_res$b_mat[,offset_var]
    ))
    mat <- mat/library_mat
  }

  if(all(is.null(gene_names))){
    if(only_nonzero) {
      idx <- which(mat != 0)
    } else {
      idx <- 1:prod(dim(mat))
    }
  } else {
    n <- nrow(mat)
    gene_idx <- which(colnames(mat) %in% gene_names)
    idx <- unlist(lapply(gene_idx, function(j){(j-1)*n+c(1:n)}))
  }

  if(length(idx) > max_num){
    idx <- sample(idx, size = max_num)
  }

  tmp_mat <- cbind(mat[idx], mean_mat[idx])
  if(bool_logscale){
    tmp_mat <- log1p(tmp_mat)
  }
  x_vec <- tmp_mat[,2]
  y_vec <- tmp_mat[,1]
  if(asp){
    if(all(is.na(xlim))) xlim <- range(c(0, x_vec, y_vec))
    if(all(is.na(ylim))) ylim <- range(c(0, x_vec, y_vec))
  } else {
    if(all(is.na(xlim))) xlim <- range(c(0, x_vec))
    if(all(is.na(ylim))) ylim <- range(c(0, y_vec))
  }

  bool_vec <- apply(tmp_mat, 1, function(x){all(x >= xlim[1]) & all(x <= xlim[2])})
  tmp_mat <- tmp_mat[which(bool_vec),]
  angle_val <- .compute_principal_angle(tmp_mat)

  graphics::plot(NA,
                 xlim = xlim,
                 ylim = ylim,
                 asp = asp,
                 ...)

  # plot diagonal
  seq_vec <- c(0, 2*max(abs(xlim)))
  graphics::points(x = x_vec,
                   y = y_vec,
                   pch = 16,
                   cex = cex,
                   col = col_points)
  graphics::lines(x = seq_vec,
                  y = seq_vec,
                  col = "black",
                  lwd = 2,
                  lty = 2)
  graphics::lines(x = seq_vec,
                  y = seq_vec * tan(angle_val*pi/180),
                  col = "green",
                  lwd = 2,
                  lty = 2)
  invisible()
}


#############
#' Compute principal angle
#'
#' @param tmp_mat a matrix with \code{n} rows (for \code{n} samples) and \code{2} columns,
#' where the first column represents the observed data and the second column represents its
#' corresponding predicted values
#'
#' @return numeric
.compute_principal_angle <- function(tmp_mat){
  pca_res <- stats::prcomp(tmp_mat, center = F, scale = F)
  vec <- pca_res$rotation[,1]; vec <- vec/eSVD2:::.l2norm(vec)
  if(sign(vec[1]) < 0)  vec <- -1*vec
  angle_val <- as.numeric(acos(as.numeric(c(0,1) %*% vec)))
  angle_val * 180/pi
}


