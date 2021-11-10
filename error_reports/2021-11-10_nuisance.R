load("2021-11-10_nuisance.RData") # load vec_list

# contains 3 different genes with the following elements:
# gene_idx: the gene index in the pseudoreal data generation
# obs_vec: the observation vector (of length equal to number of cells)
# true_vec: the true mean vector used to generate the observation vector
# est_vec: the estimated vector based on eSVD (fit on many cells and genes, but we're reporting only specific gene's vectors)
# true_nuisance: true dispersion parameter used to generate the data (i.e., r where the variance of a NB is mu+mu^2/r)
# est_nuisance: the estimated dispersion parameter using glmGamPoi::overdispersion_mle

for(i in 1:length(vec_list)){
  val2 <- glmGamPoi::overdispersion_mle(y = vec_list[[i]]$obs_vec,
                                        mean = vec_list[[i]]$est_vec)$estimate

  max_val <- max(vec_list[[i]]$true_vec)
  par(mfrow = c(1,2))
  plot(vec_list[[i]]$true_vec, vec_list[[i]]$obs_vec, pch = 16, asp = T,
       col = rgb(0.5,0.5,0.5,0.5), xlim = c(0, max_val), ylim = c(0, max_val),
       main = paste0("Gene ", vec_list[[i]]$gene_idx, " (Truth)\nNuisance: ", round(vec_list[[i]]$true_nuisance, 4)),
       xlab = "True mean", ylab = "Observed value")
  lines(c(-1e5, 1e5), c(-1e5, 1e5), col = 2, lty = 2)

  plot(vec_list[[i]]$est_vec, vec_list[[i]]$obs_vec, pch = 16, asp = T,
       col = rgb(0.5,0.5,0.5,0.5), xlim = c(0, max_val), ylim = c(0, max_val),
       main = paste0("Gene ", vec_list[[i]]$gene_idx, " (Estimated)\nNuisance: ", round(val2, 4)),
       xlab = "Estimated mean", ylab = "Observed value")
  lines(c(-1e5, 1e5), c(-1e5, 1e5), col = 2, lty = 2)
}

#########################

# ways to estimate dispersion
i <- 1
obs_vec <- vec_list[[i]]$obs_vec
est_vec <- vec_list[[i]]$est_vec

glmGamPoi::overdispersion_mle(y = obs_vec,
                              mean = est_vec)$estimate
MASS::theta.ml(y = obs_vec, mu = est_vec)
MASS::theta.mm(y = obs_vec, mu = est_vec, dfr = length(obs_vec)-1)
