rm(list=ls())
load("../../../../out/writeup8d/writeup8d_sns_pseudoreal.RData")

j <- which.max(esvd_res2$nuisance_param_vec)
esvd_res2$nuisance_param_vec[j]
true_esvd$nuisance_param_vec[j]

obs_vec <- mat[,j]
quantile(obs_vec[obs_vec>0])
length(which(obs_vec == 0))/length(obs_vec)

true_vec <- s_vec*lambda_mat[,j]
quantile(true_vec[obs_vec>0])
quantile(true_vec[obs_vec==0])
cor(obs_vec, true_vec)

nat_vec1 <- esvd_res2$x_mat %*% esvd_res2$y_mat[j,]
nat_vec2 <- esvd_res2$covariates %*% esvd_res2$b_mat[j,]
est_vec <- as.numeric(exp(nat_vec1 + nat_vec2))
quantile(est_vec[obs_vec>0])
quantile(est_vec[obs_vec==0])
cor(obs_vec, est_vec)
cor(obs_vec[obs_vec>0], est_vec[obs_vec>0])

glmgampoi_res <- glmGamPoi::overdispersion_mle(y = obs_vec,
                                               mean = est_vec)
glmgampoi_res$estimate


mean(obs_vec)
mean(est_vec)
var(obs_vec)
var(est_vec)
mean(est_vec)+mean(est_vec)^2/glmgampoi_res$estimate

png("../../../../out/fig/writeup8d/sns_pseudoreal_gene_vec.png",
    height = 1200, width = 2000, units = "px", res = 300)
set.seed(10)
par(mfrow = c(1,2))
plot(true_vec, obs_vec, pch = 16, asp = T,
     main = paste0("Gene ", j, " (Truth)"))
lines(c(-1e5, 1e5), c(-1e5, 1e5), col = 2, lty = 2)

plot(est_vec, obs_vec, pch = 16, asp = T,
     main = paste0("Gene ", j, " (Estimated)"))
lines(c(-1e5, 1e5), c(-1e5, 1e5), col = 2, lty = 2)
graphics.off()

