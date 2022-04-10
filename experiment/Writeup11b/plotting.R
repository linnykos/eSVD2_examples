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
