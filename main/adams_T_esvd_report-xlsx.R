rm(list=ls())
library(Seurat)
library(eSVD2)
library(xlsx)

load("../../../out/main/adams_T_esvd.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

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

xlsx::write.xlsx(df,
           file = "../../../out/main/adams_T_esvd_results.xlsx",
           sheetName = "adams",
           col.names = TRUE,
           row.names = TRUE,
           append = FALSE)

print("Done! :)")
