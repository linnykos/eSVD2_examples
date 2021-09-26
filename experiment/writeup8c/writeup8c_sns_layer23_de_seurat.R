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

## according to https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7678724/bin/NIHMS1053005-supplement-supplement.pdf,
# we want to regress out the following:
# [[age,  sex,  cortical  region, RIN  and  post-mortem  interval,
# as  well  as  10X  capture  and  sequencing batch and per-cell ribosomal RNA fraction]]
# a lot of these are categorical, so let's make indicators.
# We'll do a full set of indicators

categorical_var <- c("individual", "region", "sex", "Capbatch", "Seqbatch")
numerical_var <- c("age", "RNA.Integrity.Number", "post.mortem.hours", "percent.mt")
new_indicator_var <- c()
n <- ncol(sns)

for(variable in categorical_var){
  covariate <- sns@meta.data[,variable]
  uniq_level <- unique(covariate)
  for(i in uniq_level){
    tmp <- rep(0, n)
    tmp[which(covariate == i)] <- 1

    var_name <- paste0(variable, "_", i)
    sns[[var_name]] <- tmp
    new_indicator_var <- c(new_indicator_var, var_name)
  }
}

vars_to_regress <- c(numerical_var, new_indicator_var, "nFeature_RNA")
Seurat::DefaultAssay(sns) <- "RNA"
set.seed(10)
sns <- Seurat::SCTransform(sns,
                           vars.to.regress = vars_to_regress,
                           verbose = T)

set.seed(10)
sns_de <- Seurat::FindMarkers(sns,
                              assay = "SCT",
                              slot = "data",
                              ident.1 = "Control",
                              ident.2 = "ASD",
                              group.by = "diagnosis")

save(sns, sns_de,
     file = "../../../../out/writeup8d/writeup8d_sns_layer23_de_seurat.RData")

idx <- rownames(sns_de)[sns_de$p_val_adj <= 0.05]
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
de_genes <- de_genes[de_genes %in% rownames(sns_de)]
length(de_genes)
length(intersect(de_genes, idx))
sns_de[rownames(sns_de) %in% de_genes,]
