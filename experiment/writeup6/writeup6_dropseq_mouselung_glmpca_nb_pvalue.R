rm(list=ls())

load("../../../../out/writeup6/writeup6_dropseq_mouselung_glmpca_nb.RData")

stats::cor(as.numeric(glmpca_res$offsets), log(sparseMatrixStats::colSums2(mat)))

U <- as.matrix(cbind(glmpca_res$X, glmpca_res$factors))
V <- as.matrix(cbind(glmpca_res$coefX, glmpca_res$loadings))
# ignore glmpca_res$offsets since these contain the size factors

pred_mat_unnorm <- glmpca_res$glmpca_family$linkinv(tcrossprod(V,U))
pred_mat_unnorm <- t(pred_mat_unnorm)

############

membership_vec <- as.factor(lung@meta.data$celltype.l2)
celltype_vec <- names(which(table(membership_vec) > nrow(pred_mat_unnorm)/100))

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

pval_res <- pvalue_pairwise(pred_mat_unnorm, membership_vec, celltype_vec)

save.image("../../../../out/writeup6/writeup6_dropseq_mouselung_glmpca_nb_pvalue.RData")
