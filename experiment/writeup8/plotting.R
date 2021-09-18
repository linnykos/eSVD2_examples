plot_scatterplot_poisson <- function(mat,
                                     mean_mat,
                                     xlim = NA,
                                     quantile_vec = c(0.25, 0.75),
                                     only_nonzero = T,
                                     max_num = 1e5,
                                     point_cex = 1,
                                     point_col = rgb(0,0,0,0.1),
                                     mean_lwd = 2,
                                     polygon_density = 30,
                                     asp = T, ...){
  if(only_nonzero) {
    idx <- which(mat != 0)
  } else {
    idx <- 1:prod(dim(mat))
  }
  if(length(idx) > max_num){
    idx <- sample(idx, size = max_num)
  }

  if(all(is.na(xlim))) xlim <- range(c(0, mat[idx], mean_mat[idx]))
  mean_vec <- seq(0, xlim[2], length.out = 100)
  lower_vec <- stats::qpois(quantile_vec[1], lambda = mean_vec)
  upper_vec <- stats::qpois(quantile_vec[2], lambda = mean_vec)

  graphics::plot(NA,
                 xlim = xlim, ylim = xlim,
                 asp = asp,
                 ...)
  graphics::polygon(x = c(mean_vec, rev(mean_vec)),
                    y = c(lower_vec, rev(upper_vec)),
                    col = grDevices::rgb(1,0,0,0.2),
                    border = NA,
                    density = polygon_density,
                    angle = -45)
  graphics::lines(x = mean_vec, y = mean_vec,
                  col = "red",
                  lwd = mean_lwd)
  graphics::points(x = mean_mat[idx],
                   y = mat[idx],
                   pch = 16,
                   col = point_col,
                   cex = point_cex)

  invisible()
}

plot_scatterplot_nb <- function(mat,
                                prob_mat = NA,
                                mean_mat = NA,
                                size_vec,
                                main = "",
                                include_percentage_in_main = T,
                                xlim = NA,
                                quantile_vec = c(0.25, 0.75),
                                only_nonzero = T,
                                max_num = 1e5,
                                log_scale = T,
                                included_col = "black",
                                excluded_col = "red",
                                cex = 1,
                                asp = T,
                                ...){
  stopifnot(length(size_vec) == ncol(prob_mat),
            all(dim(mat) == dim(prob_mat)),
            all(is.na(mean_mat)) | all(is.na(prob_mat)))
  if(only_nonzero) {
    idx <- which(mat != 0)
  } else {
    idx <- 1:prod(dim(mat))
  }
  if(length(idx) > max_num){
    idx <- sample(idx, size = max_num)
  }

  size_mat <- matrix(size_vec,
                     nrow = nrow(mat),
                     ncol = ncol(mat),
                     byrow = T)
  if(all(is.na(mean_mat))){
    mean_mat <- eSVD2:::.mult_mat_vec(prob_mat/(1-prob_mat), size_vec)
  } else {
    prob_mat <- 1/(1+size_mat/mean_mat)
  }
  print(quantile(prob_mat))
  lower_vec <- stats::qnbinom(quantile_vec[1],
                              size = size_mat[idx],
                              prob = 1-prob_mat[idx])
  upper_vec <- stats::qnbinom(quantile_vec[2],
                              size = size_mat[idx],
                              prob = 1-prob_mat[idx])


  col_vec <- sapply(1:length(idx), function(i){
    if(mat[idx[i]] >= lower_vec[i] & mat[idx[i]] <= upper_vec[i]){
      return(included_col)
    } else {
      return(excluded_col)
    }
  })
  observed_percentage <- round(length(which(col_vec == included_col))/length(col_vec), 2)

  if(log_scale){
    x_vec <- log(mean_mat[idx])
    y_vec <- log(mat[idx])
  } else {
    x_vec <- mean_mat[idx]
    y_vec <- mat[idx]
  }

  if(include_percentage_in_main){
    main_modified = paste0(main, " (", 100*observed_percentage,  "% of ", 100*round(abs(diff(quantile_vec)),2), "%)")
  } else {
    main_modified = main
  }
  if(all(is.na(xlim))) xlim <- range(c(0, x_vec, y_vec))
  graphics::plot(x = x_vec,
                 y = y_vec,
                 xlim = xlim,
                 ylim = xlim,
                 asp = asp,
                 pch = 16,
                 cex = cex,
                 col = col_vec,
                 main = main_modified,
                 ...)
}
