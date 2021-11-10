# rm(list=ls())
# load("../../../../out/writeup8d/writeup8d_sns_pseudoreal.RData")
#
# quantile(esvd_res2$nuisance_param_vec)
# idx <- sapply(c(.5, .75, 1), function(prob){
#   val <- quantile(esvd_res2$nuisance_param_vec, prob = prob)
#   which.min(abs(esvd_res2$nuisance_param_vec - val))
# })
#
# vec_list <- lapply(idx, function(j){
#   obs_vec <- mat[,j]
#   true_vec <- s_vec*lambda_mat[,j]
#
#   nat_vec1 <- esvd_res2$x_mat %*% esvd_res2$y_mat[j,]
#   nat_vec2 <- esvd_res2$covariates %*% esvd_res2$b_mat[j,]
#   est_vec <- as.numeric(exp(nat_vec1 + nat_vec2))
#
#   true_nuisance <- true_esvd$nuisance_param_vec[j]
#   est_nuisance <- esvd_res2$nuisance_param_vec[j]
#
#   list(gene_idx = j,
#        obs_vec = obs_vec, true_vec = true_vec, est_vec = est_vec,
#        true_nuisance = true_nuisance, est_nuisance = est_nuisance)
# })
#
#
# save(vec_list,
#      file = "../../../../out/writeup8d/tmp.RData")
# j <- which.max(esvd_res2$nuisance_param_vec)
# esvd_res2$nuisance_param_vec[j]
# true_esvd$nuisance_param_vec[j]
#
#
# quantile(obs_vec[obs_vec>0])
# length(which(obs_vec == 0))/length(obs_vec)
#
#
# quantile(true_vec[obs_vec>0])
# quantile(true_vec[obs_vec==0])
# cor(obs_vec, true_vec)
#
#
#
#

#########################

load("../../out/writeup8d/tmp.RData")

# plot(log1p(est_vec), log1p(obs_vec), asp = T)
# plot(log1p(true_vec), log1p(obs_vec), asp = T)
#
# quantile(est_vec[obs_vec>0])
# quantile(est_vec[obs_vec==0])
# cor(obs_vec, est_vec)
# cor(obs_vec[obs_vec>0], est_vec[obs_vec>0])

for(i in 1:length(vec_list)){
  png(paste0("../../out/fig/writeup8d/sns_pseudoreal_gene", i, "_vec.png"),
      height = 1200, width = 2000, units = "px", res = 300)

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
  graphics.off()
}


###########################
i <- 1
obs_vec <- vec_list[[i]]$obs_vec
est_vec <- vec_list[[i]]$est_vec
true_vec <- vec_list[[i]]$true_vec

cor(obs_vec, est_vec)

glmGamPoi::overdispersion_mle(y = obs_vec,
                              mean = true_vec)$estimate
cor(obs_vec, true_vec)

idx <- which(obs_vec > 0)
glmGamPoi::overdispersion_mle(y = obs_vec[idx],
                              mean = est_vec[idx])$estimate
cor(obs_vec[idx], est_vec[idx])

glmGamPoi::overdispersion_mle(y = obs_vec[idx],
                              mean = true_vec[idx])$estimate
cor(obs_vec[idx], true_vec[idx])

idx <- which(obs_vec == 0)
glmGamPoi::overdispersion_mle(y = obs_vec[idx],
                              mean = est_vec[idx])$estimate

glmGamPoi::overdispersion_mle(y = obs_vec[idx],
                              mean = true_vec[idx])$estimate

###########################

set.seed(10)
glmGamPoi::overdispersion_mle(y = obs_vec,
                              mean = est_vec,
                              subsample = T)$estimate

set.seed(10)
glmGamPoi::overdispersion_mle(y = obs_vec,
                              mean = true_vec,
                              subsample = T)$estimate

glmGamPoi::overdispersion_mle(y = log1p(obs_vec),
                              mean = log1p(est_vec))$estimate
cor(obs_vec, est_vec)

glmGamPoi::overdispersion_mle(y = log1p(obs_vec),
                              mean = log1p(true_vec))$estimate

##############

names(obs_vec) <- NULL
names(true_vec) <- NULL
n <- length(obs_vec)
MASS::theta.ml(y = obs_vec, mu = true_vec)
MASS::theta.ml(y = obs_vec, mu = est_vec)

MASS::theta.mm(y = obs_vec, mu = true_vec, dfr = n-1)
MASS::theta.mm(y = obs_vec, mu = est_vec, dfr = n-1)

max_val <- quantile(obs_vec[obs_vec > 0], prob = 0.95)
idx <- which(obs_vec < max_val)
MASS::theta.mm(y = obs_vec[idx], mu = true_vec[idx], dfr = n-1)
MASS::theta.mm(y = obs_vec[idx], mu = est_vec[idx], dfr = n-1)

MASS::theta.mm(y = log1p(obs_vec), mu = log1p(true_vec), dfr = n-1)
MASS::theta.mm(y = log1p(obs_vec), mu = log1p(est_vec), dfr = n-1)

###################################

plot(est_vec, obs_vec, asp = T, xlim = c(0,100), ylim = c(0,100))

idx <- intersect(which(obs_vec < 100), which(obs_vec >= 1))
glmGamPoi::overdispersion_mle(y = obs_vec[idx],
                              mean = est_vec[idx])$estimate
plot(est_vec[idx], obs_vec[idx], asp = T)
plot(est_vec[idx], est_vec[idx]-obs_vec[idx], asp = T)

MASS::theta.mm(y = obs_vec[idx], mu = true_vec[idx], dfr = length(idx)-1)
MASS::theta.mm(y = obs_vec[idx], mu = est_vec[idx], dfr = length(idx)-1)


