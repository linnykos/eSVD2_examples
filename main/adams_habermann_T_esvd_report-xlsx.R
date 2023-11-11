rm(list=ls())
library(Seurat)
library(eSVD2)
library(Rmpfr)
library(openxlsx)

sheet_name_vec <- c("adams_T", "habermann_T")
filename_prefix <- "../../../out/main/"
filename_suffix <- "_esvd.RData"

wb <- openxlsx::createWorkbook()

for(sheet_name in sheet_name_vec){
  print(sheet_name)

  load(paste0(filename_prefix, sheet_name, filename_suffix))

  eSVD_obj <- eSVD2:::compute_pvalue(eSVD_obj)

  gene_names <- names(eSVD_obj$teststat_vec)
  lfc <- eSVD_obj$case_mean[gene_names] - eSVD_obj$control_mean[gene_names]
  test_stat <- eSVD_obj$teststat_vec[gene_names]
  log10pvalue <- eSVD_obj$pvalue_list$log10pvalue[gene_names]
  fdr <- eSVD_obj$pvalue_list$fdr_vec[gene_names]
  bool <- rep(FALSE, length(gene_names))
  names(bool) <- gene_names
  bool[which(fdr <= 0.05)] <- TRUE

  df <- data.frame(gene = gene_names,
                   log_fold_change = lfc,
                   test_statistic = test_stat,
                   negative_log10_pvalue = log10pvalue,
                   bh_adjusted_pvalue = fdr,
                   de_gene = bool)
  df <- df[order(df$bh_adjusted_pvalue, decreasing = F),]


  openxlsx::addWorksheet(wb = wb, sheetName = sheet_name)
  writeData(wb = wb, sheet = sheet_name, x = df)
}

openxlsx::saveWorkbook(wb,
                       file = "../../../out/main/adams-habermann-T_esvd_results.xlsx",
                       overwrite = TRUE)


