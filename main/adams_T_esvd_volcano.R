rm(list=ls())
library(Seurat)
library(eSVD2)
library(Rmpfr)
library(locfdr)

load("../../../out/main/adams_T_esvd.RData")
date_of_run

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
habermann_df_genes <- setdiff(habermann_df_genes, adams_df_genes)
de_genes <- unique(c(adams_df_genes, habermann_df_genes))
de_genes_others <- unique(c(adams_df_genes_others, habermann_df_genes_others))
cycling_genes <- c(cc.genes$s.genes, cc.genes$g2m.genes)
hk_genes <- read.csv("../../../data/housekeeping/housekeeping.txt", header = F)[,1]
de_genes_others <- setdiff(de_genes_others, de_genes)
cycling_genes <- setdiff(cycling_genes, c(de_genes_others, de_genes))
hk_genes <- setdiff(hk_genes, c(cycling_genes, de_genes_others, de_genes))

gene_names <- names(eSVD_obj$teststat_vec)
adam_idx <- which(gene_names %in% adams_df_genes)
habermann_idx <- which(gene_names %in% habermann_df_genes)
de_other_idx <- which(gene_names %in% de_genes_others)
cycling_idx <- which(gene_names %in% cycling_genes)
hk_idx <- which(gene_names %in% hk_genes)

################

logpvalue_vec <- eSVD_obj$pvalue_list$log10pvalue
selected_genes <- names(logpvalue_vec)[order(logpvalue_vec, decreasing = T)[1:length(adam_idx)]]

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
col_vec <- rep(rgb(0.6, 0.6, 0.6), p)
names(col_vec) <- gene_names
col_vec[selected_genes2] <- orange_col
col_vec[setdiff(adams_df_genes, selected_genes2)] <- purple_col

labeled_vec <- rep(FALSE, p)
names(labeled_vec) <- gene_names
labeled_vec[intersect(adams_df_genes,selected_genes)] <- TRUE

transparent_vec <- rep(TRUE, p)
names(transparent_vec) <- gene_names
transparent_vec[c(selected_genes2, adams_df_genes, hk_genes)] <- FALSE

circle_vec <- rep(FALSE, p)
names(circle_vec) <- gene_names
circle_vec[intersect(adams_df_genes, selected_genes)] <- TRUE

hk_vec <- rep(FALSE, p)
names(hk_vec) <- names(logpvalue_vec)
hk_vec[hk_genes] <- TRUE

adam_vec <- rep(FALSE, p)
names(adam_vec) <- names(logpvalue_vec)
adam_vec[adams_df_genes] <- TRUE

df <- data.frame(
  name = gene_names,
  log10pvalue = logpvalue_vec,
  lfc = log2(eSVD_obj$case_mean) - log2(eSVD_obj$control),
  adam = adam_vec,
  col = col_vec,
  circled = circle_vec,
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
plot1 <- plot1 + ggplot2::geom_point(data = subset(df, circled == TRUE & adam == TRUE), color = purple_col, size = 3, shape = 1)
plot1 <- plot1 + ggplot2::geom_point(data = subset(df, hk == TRUE), color = "white", size = 0.5)
plot1 <- plot1 + ggplot2::geom_point(data = subset(df, hk == TRUE), color = green_col, size = 0.5, alpha = .35)
plot1 <- plot1 + ggplot2::coord_cartesian(xlim = c(-max(abs(df$lfc)), max(abs(df$lfc))))
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

ggplot2::ggsave(filename = paste0("../../../out/fig/main/adams_T_volcano_ggrepel.png"),
                plot1, device = "png", width = 5*.75, height = 7*.75, units = "in")

