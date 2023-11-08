rm(list=ls())
library(Seurat)
library(eSVD2)
library(SummarizedExperiment)
library(DESeq2)


load("../../../out/main/sns_layer23_sctransform.RData")
load("../../../out/main/sns_layer23_esvd.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

hk_genes <- read.csv("../../../data/housekeeping/housekeeping.txt", header = F)[,1]
sfari_genes <- read.csv("../../../data/SFARI/SFARI-Gene_genes_09-02-2021release_01-06-2022export.csv", header = T)[,2]
deg_df <- readxl::read_xlsx(
  path = "../../../data/bulkRNA-DEG-autism/SupplementaryTable3.xlsx",
  sheet = "DEGene_Statistics"
)
deg_df <- as.data.frame(deg_df)
bulk_de_genes <- deg_df[which(deg_df[,"WholeCortex_ASD_FDR"]<=0.05),"external_gene_name"]

eSVD_obj <- eSVD2:::compute_pvalue(eSVD_obj)

gene_names <- names(eSVD_obj$case_mean)
hk_genes2 <- hk_genes[hk_genes %in% gene_names]
sfari_genes2 <- sfari_genes[sfari_genes %in% gene_names]
bulk_de_genes2 <- bulk_de_genes[bulk_de_genes %in% gene_names]

fdr_vec <- eSVD_obj$pvalue_list$fdr_vec
esvd_selected_genes <- names(fdr_vec)[which(fdr_vec <= 0.05)]
esvd_logpvalue_vec <- eSVD_obj$pvalue_list$log10pvalue
esvd_pthres <- min(esvd_logpvalue_vec[esvd_selected_genes])

sct_fdr_val <- stats::p.adjust(de_result$p_val, method = "BH")
names(sct_fdr_val) <- rownames(de_result)
sct_selected_genes <- names(sct_fdr_val)[which(sct_fdr_val <= 1e-50)]
length(sct_selected_genes)
sct_logpvalue_vec <- -log10(de_result$p_val)
names(sct_logpvalue_vec) <- rownames(de_result)
sct_pthres <- min(sct_logpvalue_vec[sct_selected_genes])

# minor adjustments
sct_logpvalue_vec <- sct_logpvalue_vec[names(esvd_logpvalue_vec)]
sct_selected_genes <- intersect(sct_selected_genes, names(esvd_logpvalue_vec))
sct_selected_genes <- names(sct_logpvalue_vec)[which(sct_logpvalue_vec >= sct_pthres)]
esvd_selected_genes <- names(esvd_logpvalue_vec)[which(esvd_logpvalue_vec >= esvd_pthres)]

#####

selected_genes <- sort(unique(c(sct_selected_genes, esvd_selected_genes)))

yellow_col <- rgb(255, 205, 114, maxColorValue = 255)
orange_col <- rgb(235, 134, 47, maxColorValue = 255)
purple_col <- rgb(122, 49, 126, maxColorValue = 255)
blue_col <- rgb(129, 139, 191, maxColorValue = 255)
green_col <- rgb(70, 177, 70, maxColorValue = 255)
red_col <- rgb(255, 90, 90, maxColorValue = 255)

p <- length(gene_names)
col_vec <- rep(rgb(0.6, 0.6, 0.6), p)
names(col_vec) <- gene_names
col_vec[setdiff(bulk_de_genes2, selected_genes)] <- blue_col
col_vec[setdiff(sfari_genes2, selected_genes)] <- purple_col
col_vec[esvd_selected_genes] <- orange_col
col_vec[sct_selected_genes] <- yellow_col
col_vec[intersect(sct_selected_genes, esvd_selected_genes)] <- red_col

labeled_vec <- rep(FALSE, p)
names(labeled_vec) <- gene_names
labeled_vec[intersect(bulk_de_genes2, selected_genes)] <- TRUE
labeled_vec[intersect(sfari_genes2, selected_genes)] <- TRUE

transparent_vec <- rep(TRUE, p)
names(transparent_vec) <- gene_names
transparent_vec[c(selected_genes, bulk_de_genes2, sfari_genes2, hk_genes2)] <- FALSE

circle_vec <- rep(FALSE, p)
names(circle_vec) <- gene_names
circle_vec[intersect(bulk_de_genes2, selected_genes)] <- TRUE
circle_vec[intersect(sfari_genes2, selected_genes)] <- TRUE

hk_vec <- rep(FALSE, p)
names(hk_vec) <- gene_names
hk_vec[hk_genes2] <- TRUE

bulk_vec <- rep(FALSE, p)
names(bulk_vec) <- gene_names
bulk_vec[bulk_de_genes2] <- TRUE

sfari_vec <- rep(FALSE, p)
names(sfari_vec) <- gene_names
sfari_vec[sfari_genes2] <- TRUE

df <- data.frame(
  name = gene_names,
  esvd_log10pvalue = esvd_logpvalue_vec,
  sct_log10pvalue = sct_logpvalue_vec,
  bulk = bulk_vec,
  col = col_vec,
  circled = circle_vec,
  hk = hk_vec,
  label = labeled_vec,
  sfari = sfari_vec,
  transparent = transparent_vec
)
order_idx <- c(which(df$transparent == T),
               intersect(which(df$transparent == F), which(!df$name %in% selected_genes)),
               intersect(which(df$transparent == F), which(df$name %in% selected_genes)))
df <- df[order_idx,]

# Create the scatterplot using ggplot2
plot1 <- ggplot2::ggplot(df, ggplot2::aes(x = sct_log10pvalue, y = esvd_log10pvalue))
plot1 <- plot1 + ggplot2::geom_hline(yintercept = esvd_pthres, color = "white", linewidth = 2)
plot1 <- plot1 + ggplot2::geom_vline(xintercept = sct_pthres, color = "white", linewidth = 2)
plot1 <- plot1 + ggplot2::geom_hline(yintercept = esvd_pthres, linetype = "dashed", color = orange_col, linewidth = 1.5)
plot1 <- plot1 + ggplot2::geom_vline(xintercept = sct_pthres, linetype = "dashed", color = yellow_col, linewidth = 1.5)
plot1 <- plot1 + ggplot2::geom_point(ggplot2::aes(color = col, alpha = ifelse(transparent, 0, 1)), show.legend = FALSE)
plot1 <- plot1 + ggplot2::scale_alpha_continuous(range = c(0.3, 1))
plot1 <- plot1 + ggplot2::scale_color_identity()
plot1 <- plot1 + ggplot2::geom_point(data = subset(df, circled == TRUE & sfari == TRUE), color = purple_col, size = 3, shape = 1)
plot1 <- plot1 + ggplot2::geom_point(data = subset(df, circled == TRUE & bulk == TRUE), color = blue_col, size = 3, shape = 1)
plot1 <- plot1 + ggplot2::geom_point(data = subset(df, hk == TRUE), color = "white", size = 0.5)
plot1 <- plot1 + ggplot2::geom_point(data = subset(df, hk == TRUE), color = green_col, size = 0.5, alpha = .35)
plot1 <- plot1 + ggplot2::labs(x = "", y = "", title = paste0("Correlation: ", round(stats::cor(esvd_logpvalue_vec, sct_logpvalue_vec), 2)))

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

ggplot2::ggsave(filename = paste0("../../../out/fig/main/sns_layer23_volcano_cross-comparison_sctransform_ggrepel.png"),
                plot1, device = "png", width = 5*.75, height = 7*.75, units = "in")

