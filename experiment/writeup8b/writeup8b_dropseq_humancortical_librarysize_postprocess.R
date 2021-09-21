rm(list=ls())
load("../../../../out/writeup8b/writeup8b_dropseq_humancortical_librarysize.RData")

plotting_func <- function(smoothing_obj,
                          xmax,
                          filename){

  median_list <- smoothing_obj$median_list
  xlim <- range(median_list[[1]][,"x"])
  xlim[2] <- min(xlim[2], xmax)
  # ylim <- range(unlist(lapply(median_list, function(x){x[,-1]})))
  png(filename,
      height = 2000, width = 3000, units = "px", res = 300)
  par(mfrow = c(2,3))
  for(i in 1:6){
    plot(NA,
         xlim = xlim,
         ylim = range(median_list[[i]][,-1]),
         xlab = "Total cell UMI count",
         ylab = "Gene UMI count",
         main = paste0("Group ", i, " (", kmean_res$size[i],
                       " genes,\nLog-exp of ",
                       round(kmean_res$centers[i,1], 2), ")"))
    polygon(x = c(median_list[[i]][,"x"], rev(median_list[[i]][,"x"])),
            y = c(median_list[[i]][,"10%"], rev(median_list[[i]][,"90%"])),
            col = grDevices::rgb(0,0,0,0.2),
            density = 30,
            angle = -45)
    lines(x = median_list[[i]][,"x"],
          y = median_list[[i]][,"50%"],
          col = "red",
          lwd = 2)
  }
  graphics.off()

  invisible()
}

zz <- rowSums(mat)
xmax <- quantile(zz, probs = 0.95)
plotting_func(smoothing_original,
              xmax,
              "../../../../out/fig/writeup8b/writeup8b_dropseq_humancortical_expression_vs_umi_original.png")

zz <- rowSums(mat)
xmax <- quantile(zz, probs = 0.95)
plotting_func(smoothing_esvd,
              xmax,
              "../../../../out/fig/writeup8b/writeup8b_dropseq_humancortical_expression_vs_umi_esvd.png")

zz <- rowSums(mat)
xmax <- quantile(zz, probs = 0.95)
plotting_func(smoothing_esvd3,
              xmax,
              "../../../../out/fig/writeup8b/writeup8b_dropseq_humancortical_expression_vs_umi_esvd3.png")

zz <- rowSums(mat)
xmax <- quantile(zz, probs = 0.95)
plotting_func(smoothing_sctransform,
              xmax,
              "../../../../out/fig/writeup8b/writeup8b_dropseq_humancortical_expression_vs_umi_sctransform.png")
