rm(list=ls())
library(Seurat)
library(eSVD2)

celltype_names <- c("cyclingta", "entprog", "ta1", "ta2")

for(celltype in celltype_names){
  print(celltype)

  load(paste0("../../../out/main/regevEpi_", celltype, "_preprocessed.RData"))

  #######################
  # first create the inflamed version
  keep_vec <- rep(1, ncol(regevEpi))
  keep_vec[which(regevEpi$Sample_Health == "Non-inflamed")] <- 0
  regevEpi$keep <- keep_vec
  regevEpi2 <- subset(regevEpi, keep == 1)

  # take only half of the healthy subjects
  tab <- table(regevEpi2$Subject, regevEpi2$Subject_Disease)
  healthy_subj <- rownames(tab[tab[,"HC"] != 0,])
  set.seed(10)
  split1 <- sample(healthy_subj, size = round(length(healthy_subj)/2), replace = F)
  split2 <- setdiff(healthy_subj, split1)
  keep_vec <- rep(1, ncol(regevEpi2))
  # Non-inflamed analysis uses split1, Inflamed uses split2
  if(any(regevEpi2$Sample_Health == "Non-inflamed")){
    keep_vec[which(regevEpi2$Subject %in% split2)] <- 0
  } else {
    keep_vec[which(regevEpi2$Subject %in% split1)] <- 0
  }
  regevEpi2$keep <- keep_vec
  regevEpi2 <- subset(regevEpi2, keep == 1)

  #######################
  # second create the non-inflamed version
  keep_vec <- rep(1, ncol(regevEpi))
  keep_vec[which(regevEpi$Sample_Health == "Inflamed")] <- 0
  regevEpi$keep <- keep_vec
  regevEpi3 <- subset(regevEpi, keep == 1)

  # take only half of the healthy subjects
  tab <- table(regevEpi3$Subject, regevEpi3$Subject_Disease)
  healthy_subj <- rownames(tab[tab[,"HC"] != 0,])
  set.seed(10)
  split1 <- sample(healthy_subj, size = round(length(healthy_subj)/2), replace = F)
  split2 <- setdiff(healthy_subj, split1)
  keep_vec <- rep(1, ncol(regevEpi3))
  # Non-inflamed analysis uses split1, Inflamed uses split2
  if(any(regevEpi3$Sample_Health == "Non-inflamed")){
    keep_vec[which(regevEpi3$Subject %in% split2)] <- 0
  } else {
    keep_vec[which(regevEpi3$Subject %in% split1)] <- 0
  }
  regevEpi3$keep <- keep_vec
  regevEpi3 <- subset(regevEpi3, keep == 1)

  ##########################

  # now compute the embeddings
  set.seed(10)
  regevEpi2 <- Seurat::RunPCA(regevEpi2, verbose = F)
  regevEpi2 <- Seurat::RunUMAP(regevEpi2, dims = 1:50)

  set.seed(10)
  regevEpi3 <- Seurat::RunPCA(regevEpi3, verbose = F)
  regevEpi3 <- Seurat::RunUMAP(regevEpi3, dims = 1:50)

  ##########################

  # now plot inflamed
  seurat_obj <- regevEpi2; status <- "inflamed"
  tab_mat <- table(seurat_obj$Subject, seurat_obj$Subject_Disease)
  case_indiv <- rownames(tab_mat)[which(tab_mat[,"Colitis"] > 0)]
  num_cases <- length(case_indiv)
  case_color_palette <- grDevices::colorRampPalette(c(rgb(140, 0, 0, maxColorValue = 255),
                                                      rgb(244, 84, 84, maxColorValue = 255)))(num_cases)
  names(case_color_palette) <- case_indiv

  control_indiv <- rownames(tab_mat)[which(tab_mat[,"HC"] > 0)]
  num_controls <- length(control_indiv)
  control_color_palette <- grDevices::colorRampPalette(c(rgb(47, 60, 190, maxColorValue = 255),
                                                         rgb(27, 198, 245, maxColorValue = 255)))(num_controls)
  names(control_color_palette) <- control_indiv
  col_palette <- c(case_color_palette, control_color_palette)

  plot1 <- Seurat::DimPlot(seurat_obj, reduction = "umap",
                           group.by = "Subject",
                           cols = col_palette)
  plot1 <- plot1 + Seurat::NoLegend()
  plot1 <- plot1 + ggplot2::ggtitle("")
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../out/fig/main/regevEpi_", celltype, "_", status, "_umap_cleaned.png"),
                  plot1, device = "png", width = 4, height = 4, units = "in",
                  dpi = 300)

  # now plot non-inflamed
  seurat_obj <- regevEpi3; status <- "noninflamed"
  tab_mat <- table(seurat_obj$Subject, seurat_obj$Subject_Disease)
  case_indiv <- rownames(tab_mat)[which(tab_mat[,"Colitis"] > 0)]
  num_cases <- length(case_indiv)
  case_color_palette <- grDevices::colorRampPalette(c(rgb(140, 0, 0, maxColorValue = 255),
                                                      rgb(244, 84, 84, maxColorValue = 255)))(num_cases)
  names(case_color_palette) <- case_indiv

  control_indiv <- rownames(tab_mat)[which(tab_mat[,"HC"] > 0)]
  num_controls <- length(control_indiv)
  control_color_palette <- grDevices::colorRampPalette(c(rgb(47, 60, 190, maxColorValue = 255),
                                                         rgb(27, 198, 245, maxColorValue = 255)))(num_controls)
  names(control_color_palette) <- control_indiv
  col_palette <- c(case_color_palette, control_color_palette)

  plot1 <- Seurat::DimPlot(seurat_obj, reduction = "umap",
                           group.by = "Subject",
                           cols = col_palette)
  plot1 <- plot1 + Seurat::NoLegend()
  plot1 <- plot1 + ggplot2::ggtitle("")
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../out/fig/main/regevEpi_", celltype, "_", status, "_umap_cleaned.png"),
                  plot1, device = "png", width = 4, height = 4, units = "in",
                  dpi = 300)
}

print("Done! :)")


