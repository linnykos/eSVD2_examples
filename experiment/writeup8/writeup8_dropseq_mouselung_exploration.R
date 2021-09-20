rm(list=ls())
load("../../../../out/writeup8/writeup8_dropseq_mouselung_esvd_initglmpca2.RData")
source("plot_library_size.R")

mat3 <- log(eSVD2:::.mult_vec_mat(1/matrixStats::rowSums2(mat2), mat2)*1e4+1)
avg_expression <- matrixStats::colMeans2(mat3)
perc_expression <- apply(mat3, 2, function(x){length(which(x > 0))/length(x)})
png("../../../../out/fig/writeup8/writeup8_dropseq_mouselung_expression.png",
    height = 1500, width = 3000, units = "px", res = 300)
par(mfrow = c(1,2))
hist(avg_expression, breaks = 50, xlab = "Mean expression (Log-normalized)", col = "gray",
     main = "Mouse lung (Dropseq)")

plot(x = avg_expression, y = perc_expression,
     pch = 16,
     xlab = "Mean expression (Log-normalized)",
     ylab = "Percent expression")
graphics.off()

set.seed(10)
kmean_res <- stats::kmeans(avg_expression, centers = 6)
cbind(kmean_res$centers, kmean_res$size)
gene_grouping <- as.factor(kmean_res$cluster)

smoothing_original <- smooth_gene_vs_umi(mat2,
                                         gene_grouping)
individual_list <- smoothing_original$individual_list
median_list <- smoothing_original$median_list

xlim <- range(median_list[[1]][,"x"])
# ylim <- range(unlist(lapply(median_list, function(x){x[,-1]})))
png("../../../../out/fig/writeup8/writeup8_dropseq_mouselung_expression_vs_umi_original.png",
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
