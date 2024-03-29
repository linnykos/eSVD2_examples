rm(list=ls())
library(Seurat)
library(eSVD2)

load("../../../out/main/habermann_T_esvd.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

df_mat <- read.csv("~/project/eSVD/data/GSE136831_adams_lung/aba1983_Data_S8.txt",
                   sep = "\t")
adams_df_genes <- df_mat$gene[which(df_mat$cellType == "T")]
adams_df_genes_others <- unique(df_mat$gene[which(df_mat$cellType %in% c("B", "Macrophage", "Macrophage Alveolar", "NK"))])
df_mat <- read.csv("~/project/eSVD/data/GSE135893_habermann_lung/Table_S4_DEG_analysis/Disease_vs_Control/T_Cells_disease_vs_control_.csv",
                   sep = ",")
habermann_df_genes <- df_mat$X
file_vec <- c("Macrophages_disease_vs_control_.csv", "Monocytes_disease_vs_control_.csv",
              "B_Cells_disease_vs_control_.csv", "NK_Cells_disease_vs_control_.csv")
habermann_df_genes_others <- unique(unlist(lapply(file_vec, function(file_suffix){
  df_mat <- read.csv(paste0("~/project/eSVD/data/GSE135893_habermann_lung/Table_S4_DEG_analysis/Disease_vs_Control/", file_suffix),
                     sep = ",")
  df_mat$X
})))
adams_df_genes <- setdiff(adams_df_genes, habermann_df_genes)
de_genes <- unique(c(adams_df_genes, habermann_df_genes))
hk_genes <- read.csv("../../../data/housekeeping/housekeeping.txt", header = F)[,1]
hk_genes <- setdiff(hk_genes, de_genes)

gene_names <- names(eSVD_obj$teststat_vec)
adam_idx <- which(gene_names %in% adams_df_genes)
habermann_idx <- which(gene_names %in% habermann_df_genes)
hk_idx <- which(gene_names %in% hk_genes)

################

logpvalue_vec <- eSVD_obj$pvalue_list$log10pvalue
selected_genes <- names(logpvalue_vec)[order(logpvalue_vec, decreasing = T)[1:length(habermann_idx)]]

###############

min_pthres <- min(logpvalue_vec[selected_genes])
selected_genes2 <- names(logpvalue_vec)[logpvalue_vec >= min_pthres]
gene_names <- names(logpvalue_vec)
adams_df_genes <- intersect(adams_df_genes, gene_names)
habermann_df_genes <- intersect(habermann_df_genes, gene_names)
hk_genes <- intersect(hk_genes, gene_names)

orange_col <- rgb(235, 134, 47, maxColorValue = 255)
purple_col <- rgb(122, 49, 126, maxColorValue = 255)
green_col <- rgb(70, 177, 70, maxColorValue = 255)

p <- length(logpvalue_vec)
lfc_vec <- log2(eSVD_obj$case_mean) - log2(eSVD_obj$control)

col_vec <- rep(rgb(0.6, 0.6, 0.6), p)
names(col_vec) <- gene_names
col_vec[selected_genes2] <- orange_col
col_vec[setdiff(habermann_df_genes, selected_genes2)] <- purple_col

xval <- quantile(lfc_vec, probs = c(0.01,0.99))
xval <- max(abs(xval))

labeled_vec <- rep(FALSE, p)
names(labeled_vec) <- gene_names
labeled_vec[intersect(
  intersect(habermann_df_genes, selected_genes),
  names(lfc_vec)[which(abs(lfc_vec) <= xval)]
)] <- TRUE

transparent_vec <- rep(TRUE, p)
names(transparent_vec) <- gene_names
transparent_vec[c(selected_genes2, habermann_df_genes, hk_genes)] <- FALSE

circle_vec <- rep(FALSE, p)
names(circle_vec) <- gene_names
circle_vec[intersect(habermann_df_genes, selected_genes)] <- TRUE

hk_vec <- rep(FALSE, p)
names(hk_vec) <- names(logpvalue_vec)
hk_vec[hk_genes] <- TRUE

habermann_vec <- rep(FALSE, p)
names(habermann_vec) <- names(logpvalue_vec)
habermann_vec[habermann_df_genes] <- TRUE

df <- data.frame(
  name = gene_names,
  log10pvalue = logpvalue_vec,
  lfc = log2(eSVD_obj$case_mean) - log2(eSVD_obj$control),
  col = col_vec,
  circled = circle_vec,
  habermann = habermann_vec,
  hk = hk_vec,
  label = labeled_vec,
  transparent = transparent_vec
)
order_idx <- c(which(df$transparent == T),
               intersect(which(df$transparent == F), which(!df$name %in% selected_genes)),
               intersect(which(df$transparent == F), which(df$name %in% selected_genes)))
df <- df[order_idx,]

# Create the scatterplot using ggplot2
plot1 <- ggplot2::ggplot(df, ggplot2::aes(x = lfc, y = log10pvalue))
plot1 <- plot1 + ggplot2::geom_hline(yintercept = min_pthres, color = "white", size = 2)
plot1 <- plot1 + ggplot2::geom_hline(yintercept = min_pthres, linetype = "dashed", color = orange_col, size = 1.5)
plot1 <- plot1 + ggplot2::geom_point(ggplot2::aes(color = col, alpha = ifelse(transparent, 0, 1)), show.legend = FALSE)
plot1 <- plot1 + ggplot2::scale_alpha_continuous(range = c(0.3, 1))
plot1 <- plot1 + ggplot2::scale_color_identity()
plot1 <- plot1 + ggplot2::geom_point(data = subset(df, circled == TRUE & habermann == TRUE), color = purple_col, size = 3, shape = 1)
plot1 <- plot1 + ggplot2::geom_point(data = subset(df, hk == TRUE), color = "white", size = 0.5)
plot1 <- plot1 + ggplot2::geom_point(data = subset(df, hk == TRUE), color = green_col, size = 0.5, alpha = .35)
plot1 <- plot1 + ggplot2::coord_cartesian(xlim = c(-xval, xval))
plot1 <- plot1 + ggplot2::labs(x = "", y = "", title = "")

# Add labels for genes where "label" is TRUE using ggrepel
plot1 <- plot1 + ggrepel::geom_text_repel(
  ggplot2::aes(label = ifelse(label, name, "")),
  box.padding = 0.25,
  force_pull = 0,
  max.overlaps = 15,
  min.segment.length = 0,
  # point.padding = 0.5,
  size = 2,  # Adjust the text size here
  segment.color = "black",  # Color of the connecting line
  segment.size = 0.5  # Thickness of the connecting line
)

ggplot2::ggsave(filename = paste0("../../../out/fig/main/habermann_T_volcano_ggrepel.png"),
                plot1, device = "png", width = 3, height = 3, units = "in",
                dpi = 600)

###############################

## https://www.pathwaycommons.org/guide/primers/statistics/fishers_exact_test/
m <- length(habermann_df_genes)
n <- length(logpvalue_vec) - m
k <- length(selected_genes2)
x <- length(intersect(selected_genes2, c(habermann_df_genes)))
fisher <- sum(sapply(x:k, function(i){
  stats::dhyper(x = i, m = m, n = n, k = k, log = F)
}))
fisher
paste0("#Author: ", m, ", #Bg: ", n, ", #Select: ", k, ", #Intersect: ", x)
