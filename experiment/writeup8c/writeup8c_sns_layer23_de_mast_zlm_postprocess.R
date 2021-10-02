load("../../../../out/writeup8c/writeup8c_sns_layer23_de_mast_zlm.RData")

# from https://github.com/himelmallick/BenchmarkSingleCell/blob/master/Library/allUtilityFunctions.R
get_pval_MAST<-function(fit){
  mat<-as.matrix(fit[,3,1])
  rownames(mat)<- row.names(fit[,,1])
  colnames(mat)<- names(metadata)

  rownames<-rownames(mat)
  colnames<-colnames(mat)
  n<-dim(mat)[2]
  List <- list()
  for(i in 1:n){
    List[[i]] <- fit[,3,2+i]
  }
  Matrix = do.call(cbind, List)
  rownames(Matrix)<-rownames
  colnames(Matrix)<-colnames
  return(Matrix)
}

coef.vector <- fit[,3,1]
pvalue.vector <- get_pval_MAST(fit)[,1]
paras <- data.frame(coef = coef.vector, pval = pvalue.vector)
paras <- paras[!duplicated(rownames(paras)),]

pval_vec <- paras$pval

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
