rm(list=ls())
load("../../../../out/writeup8c/writeup8c_sns_layer23_de_seurat.RData")

# pval_vec <- sns_de$p_val
# pval_vec <- p.adjust(pval_vec, method = "BH")
pval_vec <- sns_de$p_val_adj

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
length(de_genes)
de_genes <- de_genes[de_genes %in% rownames(sns_de)]
length(de_genes)
threshold <- 1e-30
idx <- rownames(sns_de)[which(pval_vec <= threshold)]
length(intersect(de_genes, idx))
length(de_genes)
length(idx)


sns_de[which(rownames(sns_de) %in% de_genes),]
quantile(sns_de[which(rownames(sns_de) %in% de_genes),"p_val"])
quantile(sns_de[which(rownames(sns_de) %in% de_genes),"p_val_adj"])


###############

# to look at the SCTransform results
sct_obj <- sns[["SCT"]]@SCTModel.list[[1]]
head(sct_obj@feature.attributes)

############################

sns <- Seurat::RunPCA(sns, verbose = F)
set.seed(10)
sns <- Seurat::RunUMAP(sns,
                       dims = 1:50,
                       reduction.name = 'umap.layer23',
                       reduction.key = 'rnalayerUMAP_')

covariates <- c("diagnosis", "sex", "individual", "region", "Capbatch", "Seqbatch")
for(covariate in covariates){
  plot1 <- Seurat::DimPlot(sns, reduction = "umap.layer23",
                           group.by = covariate, label = TRUE,
                           repel = TRUE, label.size = 2.5)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("SNS (Layer 2/3): ", covariate))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../../out/fig/writeup8c/sns_layer23_sctransform_umap_", covariate, ".png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
}

covariates <- c("nCount_RNA", "age", "post.mortem.hours")
for(covariate in covariates){
  plot1 <- Seurat::FeaturePlot(sns,
                               features = covariate,
                               reduction = "umap.layer23")
  plot1 <- plot1 + ggplot2::ggtitle(paste0("SNS (Layer 2/3): ", covariate))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../../out/fig/writeup8c/sns_layer23_sctransform_umap_", covariate, ".png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
}

plot1 <- Seurat::FeaturePlot(sns,
                             slot = "scale.data",
                             features = c("TTF2", "MX2", "ASCC1",
                                          "GLRA3", "CIRBP", "SAT2",
                                          "QTRT1", "CDH2", "LUC7L"),
                             reduction = "umap.layer23")
ggplot2::ggsave(filename = paste0("../../../../out/fig/writeup8c/sns_layer23_sctransform_umap_genes.png"),
                plot1, device = "png", width = 11, height = 8, units = "in")


