rm(list=ls())
library(Seurat)

load("../../../../out/writeup8d/writeup8d_sns_layer23_esvd.RData")

quantile(esvd_res2$nuisance_param_vec)

keep <- rep(0, nrow(sns@meta.data))
keep[which(rownames(sns@meta.data) %in% rownames(mat))] <- 1
sns[["keep"]] <- keep
sns <- subset(sns, keep == 1)
stopifnot(all(rownames(sns@meta.data) == rownames(mat)))

set.seed(10)
tmp <- Seurat::RunUMAP(as.matrix(esvd_res2$x_mat))@cell.embeddings
rownames(tmp) <- rownames(sns@meta.data)

sns[["esvdfactorumap"]] <- Seurat::CreateDimReducObject(embedding = tmp,
                                                        key = "esvdfactorumap_",
                                                        assay = "RNA")

covariates <- c("diagnosis", "sex", "individual", "region", "Capbatch", "Seqbatch")
for(covariate in covariates){
  plot1 <- Seurat::DimPlot(sns, reduction = "esvdfactorumap",
                           group.by = covariate, label = TRUE,
                           repel = TRUE, label.size = 2.5)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("SNS (Layer 2/3) via eSVD: ", covariate))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../../out/fig/writeup8d/sns_layer23_esvd2_umap_", covariate, ".png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
}

covariates <- c("nCount_RNA", "age", "post.mortem.hours")
for(covariate in covariates){
  plot1 <- Seurat::FeaturePlot(sns,
                               features = covariate,
                               reduction = "esvdfactorumap")
  plot1 <- plot1 + ggplot2::ggtitle(paste0("SNS (Layer 2/3) via eSVD: ", covariate))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../../out/fig/writeup8d/sns_layer23_esvd2_umap_", covariate, ".png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
}

######################

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
de_genes <- de_genes[de_genes %in% colnames(mat)]
idx <- which(colnames(mat) %in% de_genes)
colnames(esvd_res2$b_mat) <- colnames(esvd_res2$covariates)
vec <- esvd_res2$b_mat[,"diagnosis_ASD"]
xlim <- c(-2,2)
vec <- vec[vec >= xlim[1]]
vec <- vec[vec <= xlim[2]]
png("../../../../out/fig/writeup8d/sns_layer23_esvd2_hist_diagnosis.png",
    width = 1500, height = 1500, units = "px", res = 300)
hist(vec, main = "SNS (Layer 2/3) via eSVD,\nCoefficient for diagnosis (ASD)",
     col = "gray", xlab = "ASD coefficient", breaks = 50,
     xlim = xlim)
rug(vec[idx], col = "red")
graphics.off()

new_membership_vec <- rep(0, n)
asd_id <- unique(sns@meta.data[sns@meta.data$diagnosis == "ASD","individual"])
control_id <- unique(sns@meta.data[sns@meta.data$diagnosis == "Control","individual"])
for(i in 1:length(asd_id)){
  new_membership_vec[sns@meta.data$individual == asd_id[i]] <- i
}
control_id <- unique(sns@meta.data[sns@meta.data$diagnosis == "Control","individual"])
for(i in 1:length(control_id)){
  new_membership_vec[sns@meta.data$individual == control_id[i]] <- i+length(asd_id)
}

png("../../../../out/fig/writeup8d/sns_layer23_esvd2_heatmap.png",
    height = 2000, width = 2000, units = "px", res = 300)
eSVD2:::plot_scores_heatmap(esvd_res2$x_mat,
                    membership_vec = as.factor(new_membership_vec),
                    bool_log = T,
                    scaling_power = 1.5,
                    xlab = "Latent factor",
                    ylab = "Cells",
                    major_breakpoint = length(asd_id),
                    main = "SNS (Layer 2/3) via eSVD, Heatmap")
graphics.off()
