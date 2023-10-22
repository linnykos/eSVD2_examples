rm(list=ls())
library(Seurat)
library(eSVD2)
library(xlsx)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

celltype_names <- c("astpp", "endothelial", "insst", "invip", "layer4", "layer23", "layer56",
                    "layer56cc", "microglia", "oligo", "opc")

for(celltype in celltype_names){
  load(paste0("../../../out/main/sns_", celltype, "_esvd.RData"))

  # http://www.sthda.com/english/wiki/writing-data-from-r-to-excel-files-xls-xlsx
  gene <- names(eSVD_obj$pvalue_list$log10pvalue)
  mean_difference <- eSVD_obj$case_mean - eSVD_obj$control_mean
  teststat <- eSVD_obj$teststat_vec
  gaussian_teststat <- eSVD_obj$pvalue_list$gaussian_teststat
  log10pvalue <- eSVD_obj$pvalue_list$log10pvalue
  lfdr <- eSVD_obj$pvalue_list$fdr_vec

  df <- data.frame(gene = gene,
                   mean_difference = mean_difference,
                   teststat = teststat,
                   gaussian_teststat = gaussian_teststat,
                   log10pvalue = log10pvalue,
                   lfdr = lfdr)

  df <- df[order(df$log10pvalue, decreasing = T),]

  if(celltype == "astpp"){
    xlsx::write.xlsx(df,
                     file = "../../../out/main/sns_esvd_results.xlsx",
                     sheetName = celltype,
                     col.names = TRUE,
                     row.names = TRUE,
                     append = FALSE)
  } else {
    xlsx::write.xlsx(df,
                     file = "../../../out/main/sns_esvd_results.xlsx",
                     sheetName = celltype,
                     col.names = TRUE,
                     row.names = TRUE,
                     append = TRUE)

  }
}

print("Done! :)")
