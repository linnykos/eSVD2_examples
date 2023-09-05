rm(list=ls())
library(Seurat)
library(eSVD2)

load("../../../../out/main/sns_layer23_esvd.RData")
hk_genes <- read.csv("../../../../data/housekeeping/housekeeping.txt", header = F)[,1]

gene_names <- names(eSVD_obj$pvalue_list$fdr_vec)
hk_genes <- hk_genes[hk_genes %in% gene_names]

trials <- 10

# zz <- sapply(1:100, function(i){rnorm(10, mean = 1:10, sd = 1:10)})
# rowMeans(zz); apply(zz, 1, sd)

mean_mat <- eSVD_obj$fit_Second$posterior_mean_mat[,hk_genes]
var_mat <- eSVD_obj$fit_Second$posterior_var_mat[,hk_genes]

n <- nrow(mean_mat); p <- ncol(mean_mat)

pc_list <- lapply(1:trials, function(trial){
  print(trial)
  set.seed(10*trial)

  rand_mat <- matrix(
    stats::rnorm(n*p, mean = as.numeric(mean_mat), sd = sqrt(as.numeric(var_mat))),
    nrow = n,
    ncol = p
  )
  rand_mat <- scale(rand_mat)
  rand_svd <- eSVD2:::.svd_safe(
    mat = rand_mat,
    check_stability = T, # boolean
    K = 30, # positive integer
    mean_vec = NULL, # boolean, NULL or vector
    rescale = F, # boolean
    scale_max = NULL, # NULL or positive integer
    sd_vec = NULL
  )

  x_mat <- eSVD2:::.mult_mat_vec(rand_svd$u, sqrt(rand_svd$d))

  x_mat
})

###############

# anchor all other embeddings to the first one
pc_agg <- pc_list[[1]]
for(i in 2:trials){
  tmp <- crossprod(pc_list[[i]], pc_list[[1]])
  svd_tmp <- svd(tmp)
  tmp2 <- pc_list[[i]] %*% tcrossprod(svd_tmp$u, svd_tmp$v)
  pc_agg <- pc_agg + tmp2
}

rownames(pc_agg) <- colnames(sns)
colnames(pc_agg) <- paste0("eSVDPC_", 1:ncol(pc_agg))

set.seed(10)
sns[["esvdumap"]] <- Seurat::RunUMAP(pc_agg,
                                     verbose = F)

#############

tab_mat <- table(sns$individual, sns$diagnosis)
case_indiv <- rownames(tab_mat)[which(tab_mat[,"ASD"] > 0)]
num_cases <- length(case_indiv)
case_color_palette <- grDevices::colorRampPalette(c(rgb(140, 0, 0, maxColorValue = 255),
                                                    rgb(244, 84, 84, maxColorValue = 255)))(num_cases)
names(case_color_palette) <- case_indiv

control_indiv <- rownames(tab_mat)[which(tab_mat[,"Control"] > 0)]
num_controls <- length(control_indiv)
control_color_palette <- grDevices::colorRampPalette(c(rgb(47, 60, 190, maxColorValue = 255),
                                                       rgb(27, 198, 245, maxColorValue = 255)))(num_controls)
names(control_color_palette) <- control_indiv
col_palette <- c(case_color_palette, control_color_palette)

plot1 <- Seurat::DimPlot(sns, reduction = "esvdumap",
                         group.by = "individual",
                         cols = col_palette)
plot1 <- plot1 + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/fig/Writeup13b/sns_random-esvd-umap_layer23.png"),
                plot1, device = "png", width = 4, height = 4, units = "in",
                dpi = 300)
