rm(list=ls())
set.seed(10)

case_color_palette <- grDevices::colorRampPalette(c(rgb(140, 0, 0, maxColorValue = 255),
                                                    rgb(244, 84, 84, maxColorValue = 255)))(3)
control_color_palette <- grDevices::colorRampPalette(c(rgb(47, 60, 190, maxColorValue = 255),
                                                       rgb(27, 198, 245, maxColorValue = 255)))(3)
transparent_gray <- rgb(0.5,0.5,0.5,0.8)
two_letters <- substr(transparent_gray, start = 8, stop = 9)
case_color_trans_palette <- paste0(case_color_palette, two_letters)
control_color_trans_palette <- paste0(control_color_palette, two_letters)

n_each <- 50 # number of cells per person
s <- 6 # number of people
p <- 2 # number of genes
n <- n_each*s # total number of cells

covariate <- rep(1:s, each = n_each)
gene_mat <- matrix(NA, nrow = n, ncol = p)
rownames(gene_mat) <- 1:n # dummy
ind_idx <- lapply(sort(unique(covariate)), function(indiv){
  which(covariate == indiv)
})
for(i in 1:length(ind_idx)){
  rownames(gene_mat)[ind_idx[[i]]] <- paste0("c", i, "_", ind_idx[[i]])
}
case_indiv <- 1:floor(s/2)
control_indiv <- (floor(s/2)+1):s

# for gene 1
set.seed(10)
j <- 1
cc_diff <- 0.5
indiv_sd <- 0.1
group_sd <- 0.1
min_mean <- 0.5
max_mean <- 3
case_mean <- stats::runif(1, min = min_mean, max = max_mean)
control_mean <- case_mean + sample(c(-1,1), size = 1)*cc_diff

for(indiv in case_indiv){
  mean_val <- stats::rnorm(1, mean = case_mean, sd = group_sd)
  gene_mat[ind_idx[[indiv]],j] <- stats::rnorm(length(ind_idx[[indiv]]),
                                               mean = mean_val,
                                               sd = indiv_sd)
}

for(indiv in control_indiv){
  mean_val <- stats::rnorm(1, mean = control_mean, sd = group_sd)
  gene_mat[ind_idx[[indiv]],j] <- stats::rnorm(length(ind_idx[[indiv]]),
                                               mean = mean_val,
                                               sd = indiv_sd)
}

# for gene 2
set.seed(10)
j <- 2
indiv_sd <- 0.2
group_sd <- 0.1
min_mean <- 0.5
max_mean <- 2
mean_val <- stats::runif(1, min = min_mean, max = max_mean)

for(indiv in 1:length(ind_idx)){
  mean_indiv <- stats::rnorm(1, mean = mean_val, sd = group_sd)

  gene_mat[ind_idx[[indiv]],j] <- stats::rnorm(length(ind_idx[[indiv]]),
                                               mean = mean_indiv,
                                               sd = indiv_sd)
}

#####################

set.seed(10)
pca_res <- stats::prcomp(gene_mat)
denoised_matrix <- pca_res$x[,1] %*% t(pca_res$rotation[,1])
denoised_matrix <- sweep(denoised_matrix,
                         MARGIN = 2,
                         STATS = pca_res$center,
                         FUN = "+")

###############################################

pvalue_func <- function(mat, ind_idx, case_indiv, control_indiv){
  avg_posterior_mean_mat <- t(sapply(ind_idx, function(idx){
    colMeans(mat[idx,,drop=F])
  }))
  avg_posterior_var_mat <- t(sapply(ind_idx, function(idx){
    matrixStats::colVars(mat[idx,,drop=F])
  }))

  case_gaussian_mean <- colMeans(avg_posterior_mean_mat[case_indiv,,drop=F])
  control_gaussian_mean <- colMeans(avg_posterior_mean_mat[control_indiv,,drop=F])
  case_gaussian_var <- eSVD2:::.compute_mixture_gaussian_variance(
    avg_posterior_mean_mat = avg_posterior_mean_mat[case_indiv,,drop=F],
    avg_posterior_var_mat = avg_posterior_var_mat[case_indiv,,drop=F]
  )
  control_gaussian_var <- eSVD2:::.compute_mixture_gaussian_variance(
    avg_posterior_mean_mat = avg_posterior_mean_mat[control_indiv,,drop=F],
    avg_posterior_var_mat = avg_posterior_var_mat[control_indiv,,drop=F]
  )

  n1 <- length(case_indiv)
  n2 <- length(control_indiv)
  teststat_vec <- (case_gaussian_mean - control_gaussian_mean) /
    (sqrt(case_gaussian_var/n1 + control_gaussian_var/n2))

  df_num <- (case_gaussian_var/n1 + control_gaussian_var/n2)^2
  df_denom <- (case_gaussian_var/n1)^2/(n1-1) + (control_gaussian_var/n2)^2/(n2-1)
  df_vec <- df_num/df_denom

  p <- length(teststat_vec)
  pval_vec <- sapply(1:p, function(j){
    1-stats::pt(abs(teststat_vec[j]), df = df_vec[j])
  })

  list(df_vec = df_vec,
       pval_vec = pval_vec,
       teststat_vec = teststat_vec)
}

#############################

pvalue_func(gene_mat, ind_idx, case_indiv, control_indiv)
pvalue_func(denoised_matrix, ind_idx, case_indiv, control_indiv)


