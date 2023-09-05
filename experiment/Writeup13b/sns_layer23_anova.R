rm(list=ls())
library(Seurat)
library(eSVD2)

load("../../../../out/main/sns_layer23_esvd.RData")

score_mat <- eSVD_obj$fit_Second$x_mat
p <- ncol(score_mat)

tab_mat <- table(sns$individual, sns$diagnosis)
individual_vec <- sns$individual
individual_asd <- rownames(tab_mat)[which(tab_mat[,"ASD"] > 0)]
individual_control <- rownames(tab_mat)[which(tab_mat[,"Control"] > 0)]
for(individual in individual_asd){
  tmp <- which(individual_vec == individual)
  individual_vec[tmp] <- paste0("A_", individual)
}
for(individual in individual_control){
  tmp <- which(individual_vec == individual)
  individual_vec[tmp] <- paste0("C_", individual)
}
individual_vec <- as.factor(individual_vec)

pdf("../../../../out/fig/Writeup13b/sns_layer23_anova_violin.pdf",
    onefile = T, width = 8, height = 5)

for(k in 1:p){
  df <- data.frame(individual = individual_vec,
                   esvd = score_mat[,k])

  set.seed(10)
  anova_res <- stats::oneway.test(esvd ~ individual, data = df)
  # anova_res$p.value

  p1 <- ggplot2::ggplot(df, ggplot2::aes(x=individual, y=esvd))
  p1 <- p1 + ggplot2::geom_violin(trim=T, scale = "width")
  p1 <- p1 + ggplot2::geom_jitter(shape=16, position=ggplot2::position_jitter(0.2), alpha = 0.3, size = 0.5)
  p1 <- p1 + Seurat::NoLegend()
  p1 <- p1 + ggplot2::xlab("Individual")
  p1 <- p1 + ggplot2::ylab(paste0("eSVD latent dimension ", k))
  p1 <- p1 + ggplot2::geom_boxplot(width=0.05)
  p1 <- p1 + ggplot2::ggtitle(paste0("ANOVA -log10 p-value: ", round(-log10(anova_res$p.value), 2)))
  print(p1)
}

dev.off()

##############################

score_mat <- scale(eSVD_obj$fit_Second$x_mat[,-1])

set.seed(10)
sns[["esvdumap"]] <- Seurat::RunUMAP(score_mat,
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
ggplot2::ggsave(filename = paste0("../../../../out/fig/Writeup13b/sns_layer23_esvd-umap_cleaned_alt.png"),
                plot1, device = "png", width = 4, height = 4, units = "in",
                dpi = 300)
