smooth_gene_vs_umi <- function(mat, gene_grouping,
                               umi_vec = NA,
                               num_boot = 100,
                               deg = 2,
                               verbose = T){
  stopifnot(is.factor(gene_grouping))

  if(all(is.na(umi_vec))) {
    umi_vec <- matrixStats::rowSums2(mat)
  } else {
    stopifnot(length(umi_vec) == nrow(mat))
  }

  uniq_group <- sort(unique(gene_grouping))
  individual_list <- lapply(uniq_group, function(gene_group){
    idx <- which(gene_grouping == gene_group)

    np_list <- lapply(1:length(idx), function(i){
      if(verbose) print(paste0(i, " out of ", length(idx), " for group ", gene_group))
      y <- mat[,idx[i]]
      dat <- data.frame(y = y, x = umi_vec)
      np_fit <- npregfast::frfast(y ~ x,
                                  data = dat,
                                  p = deg,
                                  nboot = num_boot)

      data.frame(x = np_fit$x, y = np_fit$p[,1,1])
    })
  })

  stopifnot(length(unique(as.numeric(sapply(individual_list,
                                            function(x){range(sapply(x, nrow))}))))==1)

  median_list <- lapply(individual_list, function(lis){
    tmp <- sapply(lis, function(df){df[,"y"]})
    quantile_mat <- t(apply(tmp, 1, function(x){
      stats::quantile(x, probs = c(0.1, 0.5, 0.9))
    }))

    res <- cbind(lis[[1]][,"x"], quantile_mat)
    colnames(res)[1] <- "x"
    res
  })

  list(individual_list = individual_list,
       median_list = median_list)
}
