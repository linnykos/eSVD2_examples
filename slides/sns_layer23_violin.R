rm(list=ls())
library(Seurat)
library(eSVD2)

load("../../../out/main/sns_layer23_esvd.RData")
# load("../../../out/Writeup12/Writeup12_sns_layer23_esvd3.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

#################

df_vec <- eSVD2:::compute_df(input_obj = eSVD_obj)
teststat_vec <- eSVD_obj$teststat_vec
p <- length(teststat_vec)
gaussian_teststat <- sapply(1:p, function(j){
  qnorm(pt(teststat_vec[j], df = df_vec[j]))
})

locfdr_res <- locfdr::locfdr(gaussian_teststat, plot = 0)
fdr_vec <- locfdr_res$fdr
names(fdr_vec) <- names(gaussian_teststat)
null_mean <- locfdr_res$fp0["mlest", "delta"]
null_sd <- locfdr_res$fp0["mlest", "sigma"]
logpvalue_vec <- sapply(gaussian_teststat, function(x){
  if(x < null_mean) {
    Rmpfr::pnorm(x, mean = null_mean, sd = null_sd, log.p = T)
  } else {
    Rmpfr::pnorm(null_mean - (x-null_mean), mean = null_mean, sd = null_sd, log.p = T)
  }
})
logpvalue_vec <- -(logpvalue_vec/log(10) + log10(2))
logpvalue_vec <- pmin(logpvalue_vec, 15)

selected_genes <- names(fdr_vec)[which(fdr_vec <= 0.05)]

############

sfari_genes <- read.csv("../../../data/SFARI/SFARI-Gene_genes_09-02-2021release_01-06-2022export.csv", header = T)[,2]
deg_df <- readxl::read_xlsx(
  path = "../../../data/bulkRNA-DEG-autism/SupplementaryTable3.xlsx",
  sheet = "DEGene_Statistics"
)
deg_df <- as.data.frame(deg_df)
bulk_de_genes <- deg_df[which(deg_df[,"WholeCortex_ASD_FDR"]<=0.005),"external_gene_name"]

intersect(intersect(sfari_genes,bulk_de_genes),
          selected_genes)

gene <- "BRSK2"

indiv_tab <- table(sns$individual, sns$diagnosis)
case_indiv <- rownames(indiv_tab)[which(indiv_tab[,"ASD"] != 0)]
control_indiv <- rownames(indiv_tab)[which(indiv_tab[,"Control"] != 0)]
case_color_palette <- grDevices::colorRampPalette(c(rgb(140, 0, 0, maxColorValue = 255),
                                                    rgb(244, 84, 84, maxColorValue = 255)))(length(case_indiv))
control_color_palette <- grDevices::colorRampPalette(c(rgb(47, 60, 190, maxColorValue = 255),
                                                       rgb(27, 198, 245, maxColorValue = 255)))(length(control_indiv))
col_vec <- c(case_color_palette, control_color_palette)
names(col_vec) <- c(case_indiv, control_indiv)

df <- data.frame(value = as.numeric(sns[["RNA"]]@counts[gene,]),
                 individual = factor(sns$individual, levels = c(case_indiv, control_indiv)))

p1 <- ggplot2::ggplot(df, ggplot2::aes(x = individual, y = value, fill = individual)) +
  ggplot2::geom_violin(scale = "width", adjust = 1, draw_quantiles = c(0.25, 0.5, 0.75)) +
  ggplot2::scale_fill_manual(values = col_vec) +
  ggplot2::labs(title = "Observed counts for gene BRSK2 (15 Cases, 16 Controls)", x = "Individuals", y = "") +
  Seurat::NoLegend() + ggplot2::theme(axis.text.x = element_blank())
ggplot2::ggsave(filename = paste0("../../../out/fig/slides/sns_layer23_BRSK2_violin.png"),
                p1, device = "png", width = 10, height = 4, units = "in")
