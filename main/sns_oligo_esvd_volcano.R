rm(list=ls())
library(Seurat)
library(eSVD2)

load("../../../out/main/sns_oligo_esvd.RData")


set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

#################

fdr_vec <- eSVD_obj$pvalue_list$fdr_vec
selected_genes <- names(fdr_vec)[which(fdr_vec <= 0.05)]
gene_names <- names(fdr_vec)

############

hk_genes <- read.csv("../../../data/housekeeping/housekeeping.txt", header = F)[,1]
sfari_genes <- read.csv("../../../data/SFARI/SFARI-Gene_genes_09-02-2021release_01-06-2022export.csv", header = T)[,2]
deg_df <- readxl::read_xlsx(
  path = "../../../data/bulkRNA-DEG-autism/SupplementaryTable3.xlsx",
  sheet = "DEGene_Statistics"
)
deg_df <- as.data.frame(deg_df)
bulk_de_genes <- deg_df[which(deg_df[,"WholeCortex_ASD_FDR"]<=0.05),"external_gene_name"]

hk_genes <- hk_genes[hk_genes %in% gene_names]
sfari_genes <- sfari_genes[sfari_genes %in% gene_names]
bulk_de_genes <- bulk_de_genes[bulk_de_genes %in% gene_names]

###############

m <- length(bulk_de_genes)
n <- length(gene_names) - m
k <- length(selected_genes)
x <- length(intersect(selected_genes, bulk_de_genes))
paste0("#Author: ", m, ", #Bg: ", n, ", #Select: ", k, ", #Intersect: ", x,
       ", #Expected: ", round(m*(k/length(gene_names)),1) )
fisher <- sum(sapply(x:k, function(i){
  stats::dhyper(x = i, m = m, n = n, k = k, log = F)
}))
fisher

################

logpvalue_vec <- eSVD_obj$pvalue_list$log10pvalue
min_pthres <- min(logpvalue_vec[selected_genes])
selected_genes2 <- names(logpvalue_vec)[logpvalue_vec >= min_pthres]
gene_names <- names(eSVD_obj$teststat_vec)
bulk_de_genes <- intersect(bulk_de_genes, gene_names)
sfari_genes <- intersect(sfari_genes, gene_names)
hk_genes <- intersect(hk_genes, gene_names)

orange_col <- rgb(235, 134, 47, maxColorValue = 255)
purple_col <- rgb(122, 49, 126, maxColorValue = 255)
blue_col <- rgb(129, 139, 191, maxColorValue = 255)
green_col <- rgb(70, 177, 70, maxColorValue = 255)

p <- length(eSVD_obj$teststat_vec)
col_vec <- rep(rgb(0.6, 0.6, 0.6), p)
names(col_vec) <- gene_names
col_vec[selected_genes2] <- orange_col
col_vec[setdiff(bulk_de_genes, selected_genes2)] <- blue_col
col_vec[setdiff(sfari_genes, selected_genes2)] <- purple_col

labeled_vec <- rep(FALSE, p)
names(labeled_vec) <- gene_names
labeled_vec[intersect(bulk_de_genes,selected_genes)] <- TRUE
labeled_vec[intersect(sfari_genes,selected_genes)] <- TRUE

transparent_vec <- rep(TRUE, p)
names(transparent_vec) <- gene_names
transparent_vec[c(selected_genes2, bulk_de_genes, sfari_genes, hk_genes)] <- FALSE

circle_vec <- rep(FALSE, p)
names(circle_vec) <- gene_names
circle_vec[intersect(bulk_de_genes,selected_genes)] <- TRUE
circle_vec[intersect(sfari_genes,selected_genes)] <- TRUE

hk_vec <- rep(FALSE, p)
names(hk_vec) <- names(eSVD_obj$teststat_vec)
hk_vec[hk_genes] <- TRUE

bulk_vec <- rep(FALSE, p)
names(bulk_vec) <- names(eSVD_obj$teststat_vec)
bulk_vec[bulk_de_genes] <- TRUE

sfari_vec <- rep(FALSE, p)
names(sfari_vec) <- names(eSVD_obj$teststat_vec)
sfari_vec[sfari_genes] <- TRUE

df <- data.frame(
  name = gene_names,
  log10pvalue = logpvalue_vec,
  lfc = log2(eSVD_obj$case_mean) - log2(eSVD_obj$control),
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
plot1 <- ggplot2::ggplot(df, ggplot2::aes(x = lfc, y = log10pvalue))
plot1 <- plot1 + ggplot2::geom_hline(yintercept = min_pthres, color = "white", size = 2)
plot1 <- plot1 + ggplot2::geom_hline(yintercept = min_pthres, linetype = "dashed", color = orange_col, size = 1.5)
plot1 <- plot1 + ggplot2::geom_point(ggplot2::aes(color = col, alpha = ifelse(transparent, 0, 1)), show.legend = FALSE)
plot1 <- plot1 + ggplot2::scale_alpha_continuous(range = c(0.3, 1))
plot1 <- plot1 + ggplot2::scale_color_identity()
plot1 <- plot1 + ggplot2::geom_point(data = subset(df, circled == TRUE & sfari == TRUE), color = purple_col, size = 3, shape = 1)
plot1 <- plot1 + ggplot2::geom_point(data = subset(df, circled == TRUE & bulk == TRUE), color = blue_col, size = 3, shape = 1)
plot1 <- plot1 + ggplot2::geom_point(data = subset(df, hk == TRUE), color = "white", size = 0.5)
plot1 <- plot1 + ggplot2::geom_point(data = subset(df, hk == TRUE), color = green_col, size = 0.5, alpha = .35)
plot1 <- plot1 + ggplot2::coord_cartesian(xlim = c(-max(abs(df$lfc)), max(abs(df$lfc))))
plot1 <- plot1 + ggplot2::labs(x = "", y = "", title = paste0("Fisher -log10 p-value: ", round(-log10(fisher), 2)))

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

ggplot2::ggsave(filename = paste0("../../../out/fig/main/sns_oligo_volcano_ggrepel.png"),
                plot1, device = "png", width = 5*.75, height = 7*.75, units = "in")
