rm(list=ls())
set.seed(10)

n <- 10
p <- 13
r <- 8
k <- 3
covariate <- matrix(0, nrow = n, ncol = r)
covariate[,1] <- 1
covariate[,c(2,4)] <- rnorm(n*2)
covariate[,3] <- c(rep(1,4), rep(0,6))
covariate[,5] <- c(rep(1,4), rep(0,3), rep(1,3))
covariate[1:4,6] <- 1
covariate[5:7,7] <- 1
covariate[8:10,8] <- 1
colnames(covariate) <- c("Intercept", "Log_UMI", "Case_control", "Percent_MT",
                         "Gender", "Subject1", "Subject2", "Subject3")

coefficients <- matrix(0, nrow = p, ncol = r)
coefficients[,c(1,2,4,5)] <- rnorm(p*4)
coefficients[,3] <- c(rnorm(3), 3, 3.5, rnorm(5), -3, 0, 2.5)
coefficients[,6:8] <- rnorm(p*3, sd = 0.4)
x_mat <- matrix(rnorm(n*k), nrow = n, ncol = k)
y_mat <- matrix(rnorm(p*k), nrow = p, ncol = k)

scaling1 <- 0.5; scaling2 <- 0.5
nat_mat <- scaling1 * tcrossprod(x_mat, y_mat) + scaling2 * tcrossprod(covariate, coefficients)
obs_mat <- matrix(0, nrow = n, ncol = p)
for(i in 1:n){
  for(j in 1:p){
    obs_mat[i,j] <- rpois(1, lambda = exp(min(nat_mat[i,j], 10)))
  }
}
obs_mat <- pmin(obs_mat, 10)

###############

set.seed(10)
col_vec <- colorRampPalette(c("white",
                              rgb(98, 147, 194, maxColorValue = 255),
                              rgb(65, 88, 163, maxColorValue = 255)))(11)[-1]
tmp <- obs_mat
png("../../out/fig/main/toy_simulation_rna-heatmap.png",
    height = 2500, width = 2500, res = 500, units = "px")
par(mar = c(0.1, 0.1, 0.1, 0.1))
image(t(tmp), asp = (nrow(tmp)-1)/(ncol(tmp)-1), col = col_vec,
      xaxt = "n", yaxt = "n", bty = "n")
graphics.off()

set.seed(10)
col_vec <- colorRampPalette(c("white",
                              rgb(135, 202, 70, maxColorValue = 255)))(10)
tmp <- x_mat
png("../../out/fig/main/toy_simulation_cell_embedding.png",
    height = 2500, width = 2500, res = 500, units = "px")
par(mar = c(0.1, 0.1, 0.1, 0.1))
image(t(tmp), asp = (nrow(tmp)-1)/(ncol(tmp)-1), col = col_vec,
      xaxt = "n", yaxt = "n", bty = "n")
graphics.off()

set.seed(10)
col_vec <- colorRampPalette(c("white",
                              rgb(255, 184, 0, maxColorValue = 255)))(10)
tmp <- y_mat
png("../../out/fig/main/toy_simulation_gene_embedding.png",
    height = 2500, width = 2500, res = 500, units = "px")
par(mar = c(0.1, 0.1, 0.1, 0.1))
image(t(tmp), asp = (nrow(tmp)-1)/(ncol(tmp)-1), col = col_vec,
      xaxt = "n", yaxt = "n", bty = "n")
graphics.off()

set.seed(10)
col_vec <- colorRampPalette(c(rgb(193, 233, 154, maxColorValue = 255),
                              rgb(2, 167, 70, maxColorValue = 255)))(10)
tmp <- abs(covariate)
png("../../out/fig/main/toy_simulation_covariates.png",
    height = 2500, width = 2500, res = 500, units = "px")
par(mar = c(0.1, 0.1, 0.1, 0.1))
image(t(tmp), asp = (nrow(tmp)-1)/(ncol(tmp)-1), col = col_vec,
      xaxt = "n", yaxt = "n", bty = "n")
graphics.off()

set.seed(10)
col_vec <- colorRampPalette(c(rgb(237, 214, 154, maxColorValue = 255),
                              rgb(224, 97, 13, maxColorValue = 255)))(10)
tmp <- abs(coefficients)
png("../../out/fig/main/toy_simulation_coefficients.png",
    height = 2500, width = 2500, res = 500, units = "px")
par(mar = c(0.1, 0.1, 0.1, 0.1))
image(t(tmp), asp = (nrow(tmp)-1)/(ncol(tmp)-1), col = col_vec,
      xaxt = "n", yaxt = "n", bty = "n")
graphics.off()


