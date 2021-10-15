rm(list=ls())

library(Seurat)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

load("../../../../data/sns_autism/sns_formatted.RData")
head(sns@meta.data)
keep_vec <- rep(0, ncol(sns))
keep_vec[which(sns@meta.data$celltype == "L2/3")] <- 1
sns[["keep"]] <- keep_vec
sns <- subset(sns, keep == 1)
sns <- Seurat::NormalizeData(sns)
sns <- Seurat::FindVariableFeatures(sns,
                                    selection.method = "vst",
                                    nfeatures = 5000)

de_genes <- c("TTF2",
              "MX2",
              "ASCC1",
              "GLRA3",
              "CIRBP",
              "SAT2",
              "QTRT1",
              "CDH2",
              "LUC7L",
              "TCF25",
              "SSBP2",
              "WDR60",
              "CABP1",
              "FBLN7",
              "CDC14B",
              "GPM6A",
              "IGFBP5",
              "FAM153B",
              "GUCY1A2",
              "RAB3C",
              "SSX2IP",
              "HS6ST3",
              "TENM3",
              "DACH1",
              "PLA2G4C",
              "TOX3",
              "SPAG16",
              "FAM171B",
              "GALNTL6",
              "NUMB",
              "CAPZB",
              "DDRGK1",
              "RMST",
              "SUGP2",
              "FAM49A",
              "KCNH7",
              "BRINP3",
              "GABRB1",
              "GOLGA8B",
              "OR2L13",
              "IMMP2L",
              "ARPP19",
              "VWA8",
              "RPS15",
              "DPYSL2",
              "RFX3",
              "RSRP1",
              "NFIA",
              "SNRNP70",
              "SYN2",
              "SPIN1",
              "PLPPR4",
              "SYNPR",
              "SLC22A10",
              "LINC01378",
              "RP11-577H5.5",
              "GABRG2",
              "MIR99AHG",
              "PPP3CA",
              "MIR137HG",
              "TBRG1",
              "GGT7",
              "NLGN1",
              "GNG7",
              "FZD3",
              "LRRTM3",
              "CPE",
              "KCNJ3",
              "AQP4-AS1",
              "TRAF3",
              "PKIA",
              "MGAT4C",
              "HNRNPDL",
              "SLITRK4",
              "BMPR1B",
              "AHI1",
              "CDH9",
              "RAPGEFL1",
              "RPL34P18",
              "LINC00657",
              "COL26A1",
              "CNTN3",
              "FRMD6",
              "RP11-30J20.1",
              "SLITRK5",
              "SLC39A10",
              "STX1A",
              "RPLP2",
              "MAP2",
              "CES4A",
              "NEGR1",
              "SORBS1",
              "COL24A1",
              "VSTM2L",
              "ERBB4",
              "STARD4-AS1",
              "MAPK1",
              "HSP90AA1",
              "RPL34",
              "CNTNAP2",
              "EIF1",
              "OLFM3",
              "GRID2",
              "CHL1",
              "RAP1GAP",
              "CAMK2N1",
              "SERINC1",
              "RGS12",
              "ATP1B1")
gene_keep <- unique(c(de_genes, Seurat::VariableFeatures(sns)))
idx <- which(rownames(sns[["RNA"]]@counts) %in% gene_keep)

mat <- as.matrix(Matrix::t(sns[["RNA"]]@counts[idx,]))
print(dim(mat))
# gene_names <- colnames(mat)
# gene_names <- gene_names[-grep("^MT-", gene_names)]
# mat <- mat[,gene_names]

categorical_var <- c("diagnosis", "region", "sex", "Seqbatch") #, "individual")
numerical_var <- c("age", "RNA.Integrity.Number", "post.mortem.hours", "percent.mt", "nFeature_RNA")
n <- ncol(sns)
covariates <- as.matrix(sns@meta.data[,numerical_var])
covariates <- cbind(1, log(matrixStats::rowMeans2(mat)), covariates)
colnames(covariates)[1:2] <- c("Intercept", "Log-UMI")

for(variable in categorical_var){
  vec <- sns@meta.data[,variable]
  uniq_level <- unique(vec)
  for(i in uniq_level[-1]){
    tmp <- rep(0, n)
    tmp[which(vec == i)] <- 1

    var_name <- paste0(variable, "_", i)
    covariates <- cbind(covariates, tmp)
    colnames(covariates)[ncol(covariates)] <- var_name
  }
}

# cols_regress_out <- grep("individual", colnames(covariates))
# covariates_new <- covariates[,-cols_regress_out,drop = F]
# for(i in 1:length(cols_regress_out)){
#   df_tmp <- data.frame(covariates[,cols_regress_out[i]], covariates_new)
#   colnames(df_tmp)[1] <- "tmp"
#   lm_fit <- stats::lm("tmp ~ . - 1", data = df_tmp)
#   vec_tmp <- stats::residuals(lm_fit)
#   if(sum(abs(vec_tmp)) < 1e-6) break()
#   covariates_new <- cbind(covariates_new, vec_tmp)
#   colnames(covariates_new)[ncol(covariates_new)] <- paste0("individual_",i)
# }
#
# K <- 10
# covariates <- covariates_new[,1:(ncol(covariates_new)-2*K)]

##############
# initialization
K <- 10
n <- nrow(mat)
p <- ncol(mat)


time_start1 <- Sys.time()
init_res <- eSVD2::initialize_esvd(mat,
                                   k = K,
                                   family = "neg_binom2",
                                   covariates = covariates,
                                   column_set_to_one = "Log-UMI",
                                   offset_vec = rep(0, nrow(mat)),
                                   verbose = 1)
time_end1 <- Sys.time()
save.image("../../../../out/writeup8c/writeup8c_sns_layer23_esvd2.RData")

###################3

print("Estimating NB via eSVD, round 1")
time_start2 <- Sys.time()
set.seed(10)
esvd_res <- eSVD2::opt_esvd(init_res$x_mat,
                            init_res$y_mat,
                            mat,
                            family = "neg_binom2",
                            nuisance_param_vec = init_res$nuisance_param_vec,
                            library_size_vec = 1,
                            method = "newton",
                            b_init = init_res$b_mat,
                            covariates = init_res$covariates,
                            offset_vec = init_res$offset_vec,
                            reestimate_nuisance = T,
                            global_estimate = T,
                            reparameterize = T,
                            max_iter = 50,
                            verbose = 1)
time_end2 <- Sys.time()
save.image("../../../../out/writeup8c/writeup8c_sns_layer23_esvd2.RData")

print("Estimating NB via eSVD, round 2")
time_start3 <- Sys.time()
set.seed(10)
esvd_res2 <- eSVD2::opt_esvd(esvd_res$x_mat,
                             esvd_res$y_mat, mat,
                             family = "neg_binom2",
                             nuisance_param_vec = esvd_res$nuisance_param_vec,
                             library_size_vec = 1,
                             method = "newton",
                             b_init = esvd_res$b_mat,
                             covariates = esvd_res$covariates,
                             offset_vec = esvd_res$offset_vec,
                             reestimate_nuisance = T,
                             global_estimate = F,
                             reparameterize = T,
                             max_iter = 50,
                             tol = 1e-8,
                             verbose = 1)
time_end3 <- Sys.time()
save.image("../../../../out/writeup8c/writeup8c_sns_layer23_esvd2.RData")

