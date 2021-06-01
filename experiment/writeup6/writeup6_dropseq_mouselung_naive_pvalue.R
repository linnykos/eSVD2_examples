# from https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis.html
rm(list=ls())

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()
print("Loading in data")
dat <- anndata::read_h5ad("../../../../data/dropseq_mouselung/lung_regeneration_after_bleo")

print("Starting Seurat")
lung <- Seurat::CreateSeuratObject(counts = Matrix::t(dat$X))
lung[["celltype"]] <- dat$obs$clusters
lung <- Seurat::NormalizeData(lung, normalization.method = "LogNormalize", scale.factor = 10000)
lung <-  Seurat::FindVariableFeatures(lung, selection.method = "vst", nfeatures = 2000)

#######

mat <- t(as.matrix(lung[["RNA"]]@data[Seurat::VariableFeatures(lung),]))
membership_vec <- as.factor(lung@meta.data$celltype)
celltype_vec <- names(which(table(membership_vec) > nrow(mat)/100))

############

pvalue_pairwise <- function(mat, membership_vec, celltype_vec,
                            verbose = T){
  stopifnot(is.factor(membership_vec), all(celltype_vec %in% levels(membership_vec)))

  k <- length(celltype_vec)
  list_idx <- lapply(celltype_vec, function(x){
    which(membership_vec == x)
  })

  permn_mat <- gtools::permutations(k, 2)
  res <- matrix(NA, nrow = nrow(permn_mat), ncol = ncol(mat))
  colnames(res) <- colnames(mat)
  rownames(res) <- apply(permn_mat, 1, function(x){paste0(x, collapse = "-")})

  lookup <- data.frame(idx1 = permn_mat[,1], idx2 = permn_mat[,2],
                       type1 = celltype_vec[permn_mat[,1]],
                       type2 = celltype_vec[permn_mat[,2]])

  for(j in 1:ncol(mat)){
    print(j)
    if(verbose && ncol(mat) > 10 && j %% floor(ncol(mat)/10) == 0) cat('*')
    for(i in 1:nrow(permn_mat)){
      res[i,j] <- stats::wilcox.test(x = mat[list_idx[[permn_mat[i,1]]],j],
                                     y = mat[list_idx[[permn_mat[i,2]]],j],
                                     alternative = "greater")$p.value
    }
  }

  list(pval_mat = res, lookup_mat = lookup)
}

pval_res <- pvalue_pairwise(mat, membership_vec, celltype_vec)
save.image("../../../../out/writeup6/writeup6_dropseq_mouselung_naive_pvalue.RData")



