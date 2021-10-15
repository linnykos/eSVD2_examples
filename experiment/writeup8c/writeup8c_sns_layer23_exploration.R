rm(list=ls())

library(Seurat)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

load("../../../../data/sns_autism/sns_formatted.RData")
head(sns@meta.data)
keep_vec <- rep(0, ncol(sns))
keep_vec[which(sns@meta.data$celltype == "L2/3")] <- 1
sns[["keep"]] <- keep_vec
sns <- subset(sns, keep == 1)

set.seed(10)
sns <- Seurat::SCTransform(sns, verbose = T)
sns <- Seurat::RunPCA(sns, verbose = F)
set.seed(10)
sns <- Seurat::RunUMAP(sns, dims = 1:50,
                       reduction.name = 'umap.rna',
                       reduction.key = 'rnaUMAP_')

################

covariates <- c("diagnosis", "sex", "individual", "region", "Capbatch", "Seqbatch")
for(covariate in covariates){
  plot1 <- Seurat::DimPlot(sns, reduction = "umap.rna",
                           group.by = covariate, label = TRUE,
                           repel = TRUE, label.size = 2.5)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("SNS (Layer 2/3): ", covariate))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../../out/fig/writeup8c/sns_default_umap_", covariate, ".png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
}

covariates <- c("nCount_RNA", "age", "post.mortem.hours")
for(covariate in covariates){
  plot1 <- Seurat::FeaturePlot(sns,
                               features = covariate,
                               reduction = "umap.rna")
  plot1 <- plot1 + ggplot2::ggtitle(paste0("SNS (Layer 2/3): ", covariate))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../../out/fig/writeup8c/sns_default_umap_", covariate, ".png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
}

plot1 <- Seurat::FeaturePlot(sns,
                             features = c("TTF2", "MX2", "ASCC1",
                                          "GLRA3", "CIRBP", "SAT2",
                                          "QTRT1", "CDH2", "LUC7L"),
                             reduction = "umap.rna")
ggplot2::ggsave(filename = paste0("../../../../out/fig/writeup8c/sns_default_umap_genes.png"),
                plot1, device = "png", width = 11, height = 8, units = "in")

###################

set.seed(10)
sns_de <- Seurat::FindMarkers(sns,
                              assay = "SCT",
                              slot = "data",
                              ident.1 = "Control",
                              ident.2 = "ASD",
                              group.by = "diagnosis",
                              logfc.threshold = 0,
                              min.pct = 0,
                              min.cells.feature = 0,
                              verbose = T)

save(sns, sns_de,
     file = "../../../../out/writeup8c/writeup8c_sns_layer23_exploration.RData")

######################

sign_vec <- rep(1, nrow(sns_de))
sign_vec[which(sns_de$avg_log2FC > 0)] <- -1 #if the log-fold change is positive, control is larger than case
z_val <- stats::qnorm(0.5+(1-sns_de$p_val)/2)
z_val[is.infinite(z_val)] <- 10
z_val <- z_val * sign_vec

de_genes <- c("TTF2",
              "MX2",
              "ASCC1",
              "GLRA3",
              "CIRBP",
              "SAT2",
              "QTRT1",
              "CDH2",
              "LUC7L",
              "TCF25",
              "SSBP2",
              "WDR60",
              "CABP1",
              "FBLN7",
              "CDC14B",
              "GPM6A",
              "IGFBP5",
              "FAM153B",
              "GUCY1A2",
              "RAB3C",
              "SSX2IP",
              "HS6ST3",
              "TENM3",
              "DACH1",
              "PLA2G4C",
              "TOX3",
              "SPAG16",
              "FAM171B",
              "GALNTL6",
              "NUMB",
              "CAPZB",
              "DDRGK1",
              "RMST",
              "SUGP2",
              "FAM49A",
              "KCNH7",
              "BRINP3",
              "GABRB1",
              "GOLGA8B",
              "OR2L13",
              "IMMP2L",
              "ARPP19",
              "VWA8",
              "RPS15",
              "DPYSL2",
              "RFX3",
              "RSRP1",
              "NFIA",
              "SNRNP70",
              "SYN2",
              "SPIN1",
              "PLPPR4",
              "SYNPR",
              "SLC22A10",
              "LINC01378",
              "RP11-577H5.5",
              "GABRG2",
              "MIR99AHG",
              "PPP3CA",
              "MIR137HG",
              "TBRG1",
              "GGT7",
              "NLGN1",
              "GNG7",
              "FZD3",
              "LRRTM3",
              "CPE",
              "KCNJ3",
              "AQP4-AS1",
              "TRAF3",
              "PKIA",
              "MGAT4C",
              "HNRNPDL",
              "SLITRK4",
              "BMPR1B",
              "AHI1",
              "CDH9",
              "RAPGEFL1",
              "RPL34P18",
              "LINC00657",
              "COL26A1",
              "CNTN3",
              "FRMD6",
              "RP11-30J20.1",
              "SLITRK5",
              "SLC39A10",
              "STX1A",
              "RPLP2",
              "MAP2",
              "CES4A",
              "NEGR1",
              "SORBS1",
              "COL24A1",
              "VSTM2L",
              "ERBB4",
              "STARD4-AS1",
              "MAPK1",
              "HSP90AA1",
              "RPL34",
              "CNTNAP2",
              "EIF1",
              "OLFM3",
              "GRID2",
              "CHL1",
              "RAP1GAP",
              "CAMK2N1",
              "SERINC1",
              "RGS12",
              "ATP1B1")

idx <- which(rownames(sns_de) %in% de_genes)
png("../../../../out/fig/writeup8c/sns_default_hist_zvalues.png",
    width = 1500, height = 1500, units = "px", res = 300)
hist(z_val, main = "SNS (Layer 2/3),\nZ-value from vanilla analysis",
     col = "gray", xlab = "Z-value", breaks = 50)
rug(jitter(z_val[idx]), col = "red")
graphics.off()

# save(sns_de, file = "../../../../out/writeup8c/sns_layer23_exploration_de.RData")
# load("../../out/writeup8c/sns_layer23_exploration_de.RData")

# see https://bioconductor.org/packages/release/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html
sns_de$avg_log2FC <- -sns_de$avg_log2FC
png("../../../../out/fig/writeup8c/sns_default_volcano.png",
    width = 3000, height = 3000, units = "px", res = 300)
EnhancedVolcano:: EnhancedVolcano(sns_de,
                                  lab = rownames(sns_de),
                                  x = 'avg_log2FC',
                                  y = 'p_val',
                                  xlim = c(-1,1),
                                  selectLab = rownames(sns_de)[rownames(sns_de) %in% de_genes],
                                  FCcutoff = log2(1.10),
                                  pCutoff = 1e-30)
graphics.off()

