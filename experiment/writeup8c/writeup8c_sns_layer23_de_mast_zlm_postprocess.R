rm(list=ls())
load("../../../../out/writeup8c/writeup8c_sns_layer23_de_mast_zlm.RData")

dimnames(fit)

# see https://www.bioconductor.org/help/course-materials/2016/BioC2016/ConcurrentWorkshops2/McDavid/MAST-intro.pdf
# from https://github.com/himelmallick/BenchmarkSingleCell/blob/master/Library/allUtilityFunctions.R.
# I think this grabs the wrong values though...
# get_pval_MAST<-function(fit, sca){
#   mat <- as.matrix(fit[,3,3])
#   rownames(mat) <- row.names(fit[,,1])
#   colnames(mat) <- names(SummarizedExperiment::colData(sca))
#
#   rownames <- rownames(mat)
#   colnames <- colnames(mat)
#   n <- dim(mat)[2]
#   List <- list()
#   for(i in 1:n){
#     List[[i]] <- fit[,3,2+i]
#   }
#   Matrix = do.call(cbind, List)
#   rownames(Matrix)<-rownames
#   colnames(Matrix)<-colnames
#   return(Matrix)
# }
#
# coef.vector <- fit[,3,1]
# pvalue.vector <- get_pval_MAST(fit, sca)[,1]
# paras <- data.frame(coef = coef.vector, pval = pvalue.vector)
# paras <- paras[!duplicated(rownames(paras)),]
pval_vec <- fit[,3,3]

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
de_genes <- de_genes[de_genes %in% rownames(paras)]
length(de_genes)
threshold <- 1e-5
pval_vec2 <- p.adjust(pval_vec, method = "bonferroni")
idx <- rownames(paras)[which(pval_vec2 <= threshold)]
length(intersect(de_genes, idx))
length(de_genes)
length(idx)

#########################



# calculate the log2-fold change
# from https://github.com/satijalab/seurat/blob/master/R/differential_expression.R
mean.fxn <- function(x){log(x = rowMeans(x = expm1(x = x)) + 1, base = 2)}
fc.name <- "avg_log"
cells.1 <- rownames(sns@meta.data)[which(sns@meta.data$diagnosis == "ASD")]
cells.2 <- rownames(sns@meta.data)[which(sns@meta.data$diagnosis == "Control")]
features <- rownames(sns[["RNA"]]@counts)
object <- log1p(sns[["RNA"]]@counts)
data.1 <- mean.fxn(object[features, cells.1, drop = FALSE])
data.2 <- mean.fxn(object[features, cells.2, drop = FALSE])
fc <- (data.1 - data.2)

sns_de <- data.frame(avg_log2FC = fc, p_val = pval_vec)
rownames(sns_de) <- names(fc)

sign_vec <- rep(1, nrow(sns_de))
sign_vec[which(sns_de$avg_log2FC < 0)] <- -1
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
png("../../../../out/fig/writeup8c/sns_layer23_mast_hist_zvalues.png",
    width = 1500, height = 1500, units = "px", res = 300)
hist(z_val, main = "SNS (Layer 2/3),\nZ-value from vanilla analysis",
     col = "gray", xlab = "Z-value", breaks = 50)
rug(jitter(z_val[idx]), col = "red")
graphics.off()

# save(sns_de, file = "../../../../out/writeup8c/sns_layer23_mast.RData")
# load("../../out/writeup8c/sns_layer23_mast.RData")

png("../../out/fig/writeup8c/sns_layer23_mast_volcano.png",
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


