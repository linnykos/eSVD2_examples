rm(list=ls())

library(Seurat)
load("../../../out/main/sns_layer23_processed.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

sns$region <- factor(sns$region)
sns$sex <- factor(sns$sex)
sns$Seqbatch <- factor(sns$Seqbatch)
sns$Capbatch <- factor(sns$Capbatch)
set.seed(10)
sns <- Seurat::SCTransform(sns, method = "glmGamPoi",
                             residual.features = sns[["RNA"]]@var.features,
                             vars.to.regress = c("region", "sex", "Seqbatch", "Capbatch", "age", "RNA.Integrity.Number", "post.mortem.hours", "percent.mt"),
                             verbose = T)
Seurat::Idents(sns) <- "diagnosis"
levels(sns)

Seurat::DefaultAssay(sns) <- "SCT"
de_result <- Seurat::FindMarkers(sns, ident.1 = "ASD", ident.2 = "Control",
                                 slot = "scale.data",
                                 test.use = "wilcox",
                                 logfc.threshold = 0,
                                 min.pct = 0,
                                 verbose = T)


case_idx <- which(Seurat::Idents(sns) == "ASD")
control_idx <- which(Seurat::Idents(sns) == "Control")
n1 <- length(case_idx); n2 <- length(control_idx)
# see https://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U_test#Normal_approximation_and_tie_correction
null_mean <- n1*n2/2
null_sd <- sqrt(n1*n2*(n1+n2+1)/12)

wilcox_stats_list <- lapply(1:nrow(de_result), function(i){
  if(i %% floor(nrow(de_result)/10) == 0) cat('*')
  gene_name <- rownames(de_result)[i]

  x <- sns[["SCT"]]@scale.data[gene_name, case_idx]
  y <- sns[["SCT"]]@scale.data[gene_name, control_idx]

  wilcox_res <- stats::wilcox.test(x = x, y = y)
  test_stat <- wilcox_res$statistic
  z_score <- (test_stat-null_mean)/null_sd
  p_val_check <- wilcox_res$p.value

  c(teststatistics = test_stat,
    null_mean = null_mean,
    null_sd = null_sd,
    z_score = z_score,
    p_val_check = p_val_check)
})

wilcox_mat <- do.call(rbind, wilcox_stats_list)
colnames(wilcox_mat) <- c("teststatistics", "null_mean", "null_sd", "z_score", "p_val_check")

de_result <- cbind(de_result, wilcox_mat)

save(sns, de_result,
     date_of_run, session_info,
     file = "../../../out/main/sns_layer23_sctransform.RData")
