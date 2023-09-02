rm(list=ls())
library(Seurat)
library(eSVD2)

load("../../../out/main/regevEpi_ta1-inflamed_esvd.RData")
eSVD_obj_inflamed <- eSVD_obj
load("../../../out/main/regevEpi_ta1-noninflamed_esvd.RData")
eSVD_obj_noninflamed <- eSVD_obj

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

sheet1 <- as.data.frame(readxl::read_xlsx("~/nzhanglab/data/SCP259_regev_colitis/NIHMS1532849-supplement-11.xlsx",
                                          sheet = "Epithelial (Non-Inflamed vs. He"))
sheet2 <- as.data.frame(readxl::read_xlsx("~/nzhanglab/data/SCP259_regev_colitis/NIHMS1532849-supplement-11.xlsx",
                                          sheet = "Epithelial (Inflamed vs. Health"))
sheet3 <- as.data.frame(readxl::read_xlsx("~/nzhanglab/data/SCP259_regev_colitis/NIHMS1532849-supplement-11.xlsx",
                                          sheet = "Epithelial (Inflamed vs. Non-In"))
noninf_de_genes <- sheet1[sheet1$ident == "TA 1","gene"]
inf_de_genes <- sheet2[sheet2$ident == "TA 1","gene"]
hk_genes <- read.csv("../../../data/housekeeping/housekeeping.txt", header = F)[,1]

gene_names <- names(eSVD_obj$teststat_vec)
inf_de_idx <- which(gene_names %in% inf_de_genes)
noninf_de_idx <- which(gene_names %in% noninf_de_genes)
hk_idx <- which(gene_names %in% hk_genes)
hk_idx <- setdiff(hk_idx, c(inf_de_idx, noninf_de_idx))

orange_col <- rgb(235, 134, 47, maxColorValue = 255)
purple_col <- rgb(122, 49, 126, maxColorValue = 255)
green_col <- rgb(70, 177, 70, maxColorValue = 255)
green_col_trans <- rgb(70, 177, 70, 255*.35, maxColorValue = 255)

######################

eSVD_obj <- eSVD_obj_inflamed  ## change this line
logpvalue_vec <- eSVD_obj$pvalue_list$log10pvalue
gene_names <- names(logpvalue_vec)
noninf_de_genes <- intersect(noninf_de_genes, gene_names)
inf_de_genes <- intersect(inf_de_genes, gene_names)
author_genes <- inf_de_genes ## change this line
lfc_vec <- log2(eSVD_obj$case_mean) - log2(eSVD_obj$control)

hk_genes <- intersect(hk_genes, gene_names)
selected_genes <- names(logpvalue_vec)[order(logpvalue_vec, decreasing = T)[1:length(author_genes)]]
min_pthres <- min(logpvalue_vec[selected_genes])
selected_genes2 <- names(logpvalue_vec)[logpvalue_vec >= min_pthres]

orange_col <- rgb(235, 134, 47, maxColorValue = 255)
purple_col <- rgb(122, 49, 126, maxColorValue = 255)
green_col <- rgb(70, 177, 70, maxColorValue = 255)

p <- length(logpvalue_vec)
col_vec <- rep(rgb(0.6, 0.6, 0.6), p)
names(col_vec) <- gene_names
col_vec[selected_genes2] <- orange_col
col_vec[setdiff(author_genes, selected_genes2)] <- purple_col

xval <- max(abs(lfc_vec))

labeled_vec <- rep(FALSE, p)
names(labeled_vec) <- gene_names
labeled_vec[intersect(
  intersect(author_genes,selected_genes),
  names(lfc_vec)[which(abs(lfc_vec) <= xval)]
)] <- TRUE

transparent_vec <- rep(TRUE, p)
names(transparent_vec) <- gene_names
transparent_vec[c(selected_genes2, author_genes, hk_genes)] <- FALSE

circle_vec <- rep(FALSE, p)
names(circle_vec) <- gene_names
circle_vec[intersect(author_genes, selected_genes)] <- TRUE

hk_vec <- rep(FALSE, p)
names(hk_vec) <- names(logpvalue_vec)
hk_vec[hk_genes] <- TRUE

author_vec <- rep(FALSE, p)
names(author_vec) <- names(logpvalue_vec)
author_vec[author_genes] <- TRUE

df <- data.frame(
  name = gene_names,
  log10pvalue = logpvalue_vec,
  lfc = lfc_vec,
  col = col_vec,
  circled = circle_vec,
  author = author_vec,
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
plot1 <- plot1 + ggplot2::geom_point(data = subset(df, circled == TRUE & author == TRUE), color = purple_col, size = 3, shape = 1)
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

ggplot2::ggsave(filename = paste0("../../../out/fig/main/regevEpi_ta1-inflamed_volcano_ggrepel.png"),
                plot1, device = "png", width = 3, height = 3, units = "in",
                dpi = 600)

#####################

eSVD_obj <- eSVD_obj_noninflamed  ## change this line
logpvalue_vec <- eSVD_obj$pvalue_list$log10pvalue
gene_names <- names(logpvalue_vec)
noninf_de_genes <- intersect(noninf_de_genes, gene_names)
inf_de_genes <- intersect(inf_de_genes, gene_names)
author_genes <- noninf_de_genes ## change this line
lfc_vec <- log2(eSVD_obj$case_mean) - log2(eSVD_obj$control)

hk_genes <- intersect(hk_genes, gene_names)
selected_genes <- names(logpvalue_vec)[order(logpvalue_vec, decreasing = T)[1:length(author_genes)]]
min_pthres <- min(logpvalue_vec[selected_genes])
selected_genes2 <- names(logpvalue_vec)[logpvalue_vec >= min_pthres]

orange_col <- rgb(235, 134, 47, maxColorValue = 255)
purple_col <- rgb(122, 49, 126, maxColorValue = 255)
green_col <- rgb(70, 177, 70, maxColorValue = 255)

p <- length(logpvalue_vec)
col_vec <- rep(rgb(0.6, 0.6, 0.6), p)
names(col_vec) <- gene_names
col_vec[selected_genes2] <- orange_col
col_vec[setdiff(author_genes, selected_genes2)] <- purple_col

xval <- max(abs(lfc_vec))

labeled_vec <- rep(FALSE, p)
names(labeled_vec) <- gene_names
labeled_vec[intersect(
  intersect(author_genes,selected_genes),
  names(lfc_vec)[which(abs(lfc_vec) <= xval)]
)] <- TRUE

transparent_vec <- rep(TRUE, p)
names(transparent_vec) <- gene_names
transparent_vec[c(selected_genes2, author_genes, hk_genes)] <- FALSE

circle_vec <- rep(FALSE, p)
names(circle_vec) <- gene_names
circle_vec[intersect(author_genes, selected_genes)] <- TRUE

hk_vec <- rep(FALSE, p)
names(hk_vec) <- names(logpvalue_vec)
hk_vec[hk_genes] <- TRUE

author_vec <- rep(FALSE, p)
names(author_vec) <- names(logpvalue_vec)
author_vec[author_genes] <- TRUE

df <- data.frame(
  name = gene_names,
  log10pvalue = logpvalue_vec,
  lfc = lfc_vec,
  col = col_vec,
  circled = circle_vec,
  author = author_vec,
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
plot1 <- plot1 + ggplot2::geom_point(data = subset(df, circled == TRUE & author == TRUE), color = purple_col, size = 3, shape = 1)
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

ggplot2::ggsave(filename = paste0("../../../out/fig/main/regevEpi_ta1-noninflamed_volcano_ggrepel.png"),
                plot1, device = "png", width = 3, height = 3, units = "in",
                dpi = 600)

###################################

eSVD_obj <- eSVD_obj_inflamed  ## change this line
logpvalue_vec <- eSVD_obj$pvalue_list$log10pvalue
gene_names <- names(logpvalue_vec)
noninf_de_genes <- intersect(noninf_de_genes, gene_names)
inf_de_genes <- intersect(inf_de_genes, gene_names)
author_genes <- inf_de_genes ## change this line
selected_genes <- names(logpvalue_vec)[order(logpvalue_vec, decreasing = T)[1:length(author_genes)]]
min_pthres <- min(logpvalue_vec[selected_genes])
selected_genes2 <- names(logpvalue_vec)[logpvalue_vec >= min_pthres]

## https://www.pathwaycommons.org/guide/primers/statistics/fishers_exact_test/
m <- length(author_genes)
n <- length(logpvalue_vec) - m
k <- length(selected_genes2)
x <- length(intersect(selected_genes2, author_genes))
fisher <- sum(sapply(x:k, function(i){
  stats::dhyper(x = i, m = m, n = n, k = k, log = F)
}))
fisher
paste0("#Author: ", m, ", #Bg: ", n, ", #Select: ", k, ", #Intersect: ", x)

####

eSVD_obj <- eSVD_obj_noninflamed  ## change this line
logpvalue_vec <- eSVD_obj$pvalue_list$log10pvalue
gene_names <- names(logpvalue_vec)
noninf_de_genes <- intersect(noninf_de_genes, gene_names)
inf_de_genes <- intersect(inf_de_genes, gene_names)
author_genes <- noninf_de_genes ## change this line
selected_genes <- names(logpvalue_vec)[order(logpvalue_vec, decreasing = T)[1:length(author_genes)]]
min_pthres <- min(logpvalue_vec[selected_genes])
selected_genes2 <- names(logpvalue_vec)[logpvalue_vec >= min_pthres]

## https://www.pathwaycommons.org/guide/primers/statistics/fishers_exact_test/
m <- length(author_genes)
n <- length(logpvalue_vec) - m
k <- length(selected_genes2)
x <- length(intersect(selected_genes2, author_genes))
fisher <- sum(sapply(x:k, function(i){
  stats::dhyper(x = i, m = m, n = n, k = k, log = F)
}))
fisher
paste0("#Author: ", m, ", #Bg: ", n, ", #Select: ", k, ", #Intersect: ", x)

