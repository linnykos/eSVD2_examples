rm(list=ls())
#
# library(Seurat); library(eSVD2)
#
# set.seed(10)
# date_of_run <- Sys.time()
# session_info <- devtools::session_info()
#
# load("../../../../out/writeup7/writeup7_sns_esvd_nb_part1.RData")
#
# nat_mat <- tcrossprod(esvd_res$x_mat, esvd_res$y_mat) + tcrossprod(covariates, esvd_res$b_mat)
# print("Starting NB nuisance")
# time_start3 <- Sys.time()
# nuisance_vec <- eSVD2::initialize_nuisance_param(mat, nat_mat, family = "neg_binom",
#                                                  library_size_vec = 1)
# time_end3 <- Sys.time()
# rm(list = "nat_mat")
# save.image("../../../../out/writeup7/writeup7_sns_esvd_nb_part2.RData")
load("../../../../out/writeup7/writeup7_sns_esvd_nb_part2.RData")

print("Starting NB initialization")
time_start4 <- Sys.time()

#####################
dat <- mat
set.seed(10)
k = K
family = "neg_binom"
nuisance_param_vec = nuisance_vec
library_size_vec = 1
check_rank = F
config = eSVD2::initialization_options()
verbose = 1


stopifnot(is.character(family))
if(family != "gaussian") stopifnot(all(dat[!is.na(dat)] >= 0))
if(all(!is.na(nuisance_param_vec)) & length(nuisance_param_vec) == 1) {
  nuisance_param_vec <- rep(nuisance_param_vec[1], ncol(dat))
}

if(verbose > 0) print(paste0(Sys.time(),": Rescaling data"))
n <- nrow(dat); p <- ncol(dat)
family <- eSVD2:::.string_to_distr_funcs(family)
library_size_vec <- eSVD2:::.parse_library_size(dat, library_size_vec = library_size_vec)
rescaled_dat <- t(sapply(1:nrow(dat), function(i){
  dat[i,]/library_size_vec[i]
}))
rm(list = c("mat", "dat")); gc(T)

if(verbose > 0) print(paste0(Sys.time(),": Applying matrix completion"))
rescaled_dat <- eSVD2:::.matrix_completion(rescaled_dat, k = k)
if(verbose > 0) print(paste0(Sys.time(),": Determining initial matrix"))
init_res <- eSVD2:::.determine_initial_matrix(rescaled_dat, k = k, family = family,
                                      nuisance_param_vec = nuisance_param_vec,
                                      max_val = config$max_val,
                                      tol = config$tol)
nat_mat <- init_res$nat_mat; domain <- init_res$domain

verbose = 2
range(nat_mat)
quantile(nuisance_param_vec)
zz <- which(is.infinite(nat_mat), arr.ind = T)
