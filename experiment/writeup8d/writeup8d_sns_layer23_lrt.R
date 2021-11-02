rm(list=ls())
library(Seurat)
load("../../../../out/writeup8d/writeup8d_sns_layer23_esvd.RData")
esvd_res_full <- esvd_res2

load("../../../../out/writeup8d/writeup8d_sns_layer23_esvd_withoutdiagnosis.RData")
esvd_res_limited <- esvd_res2

ls_vec <- ls()
ls_vec <- ls_vec[!ls_vec %in% c("mat", "esvd_res_full", "esvd_res_limited")]
rm(list = ls_vec)

############

cor(esvd_res_limited$nuisance_param_vec, esvd_res_full$nuisance_param_vec)
quantile(esvd_res_limited$nuisance_param_vec)
quantile(esvd_res_full$nuisance_param_vec)

png(file = "../../../../out/fig/writeup8e/sns_layer23_lrt_nuisance.png",
    height = 1200, width = 1200, units = "px", res = 300)
plot(esvd_res_limited$nuisance_param_vec,
     esvd_res_full$nuisance_param_vec,
     pch = 16, col = rgb(0.5, 0.5, 0.5, 0.5), asp = T,
     xlab = "Limited nuisance", ylab = "Full nuisance")
graphics.off()

###########

p <- nrow(esvd_res_full$y_mat)
lrt_test_vec <- sapply(1:p, function(j){
  mean_vec_limited <- esvd_res_limited$x_mat %*% esvd_res_limited$y_mat[j,] +
                            esvd_res_limited$covariates %*% esvd_res_limited$b_mat[j,]
  mean_vec_full <- esvd_res_full$x_mat %*% esvd_res_full$y_mat[j,] +
                         esvd_res_full$covariates %*% esvd_res_full$b_mat[j,]

  r_limited <- esvd_res_limited$nuisance_param_vec[j]
  r_full <- esvd_res_full$nuisance_param_vec[j]

  ll_limited <- -sum(-(mat[,j] * mean_vec_limited) +
    mat[,j]*log(exp(mean_vec_limited) + r_limited) +
    r_limited * log(exp(mean_vec_limited) + r_limited))
  ll_full <- -sum(-(mat[,j] * mean_vec_full) +
    mat[,j]*log(exp(mean_vec_full) + r_full) +
    r_full * log(exp(mean_vec_full) + r_full))
  #[[ this seems like one weakness of fitting twice -- you're not guaranteed that the full is larger likelihood than the limited]]

  -2*(ll_limited - ll_full)

})

quantile(lrt_test_vec, probs = c(0.05, 0.95))

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
lrt_test_vec[lrt_test_vec <= -7000] <- NA
lrt_test_vec[lrt_test_vec >= 6000] <- NA
png("../../../../out/fig/writeup8e/sns_layer23_lrt_teststat.png",
    width = 1500, height = 1500, units = "px", res = 300)
hist(lrt_test_vec, main = "SNS (Layer 2/3) via eSVD,\nLRT test-stat for diagnosis (ASD)",
     col = "gray", xlab = "ASD coefficient", breaks = 50)
rug(lrt_test_vec[idx], col = "red")
graphics.off()
