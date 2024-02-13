rm(list=ls())
library(Seurat)
library(SummarizedExperiment)
library(DESeq2)
library(UniIsoRegression)

source("roc_functions.R")

for(ii in 1:3){
  print(paste0("Working on gene setting: ", ii))

  for(jj in 1:3){
    print(paste0("Working on size_factor setting: ", jj))

    load(paste0("~/kzlinlab/projects/eSVD2/out/simulation/simulation-power_geneSetting",
                ii, "_individualSetting", jj,
                ".RData"))
    load(paste0("~/kzlinlab/projects/eSVD2/out/simulation/simulation-power-esvd_geneSetting",
                ii, "_individualSetting", jj,
                ".RData"))
    load(paste0("~/kzlinlab/projects/eSVD2/out/simulation/simulation-power-deseq2_geneSetting",
                ii, "_individualSetting", jj,
                ".RData"))
    load(paste0("~/kzlinlab/projects/eSVD2/out/simulation/simulation-power-sctransform_geneSetting",
                ii, "_individualSetting", jj,
                ".RData"))
    load(paste0("~/kzlinlab/projects/eSVD2/out/simulation/simulation-power-mast_geneSetting",
                ii, "_individualSetting", jj,
                ".RData"))

    #######

    p <- length(true_fdr_vec)
    true_de_idx <- which(true_fdr_vec < 0.05)

    eSVD_obj <- eSVD2:::compute_pvalue(input_obj = eSVD_obj)
    esvd_pvalue <- 10^(-eSVD_obj$pvalue_list$log10pvalue[paste0("gene-", 1:p)])
    esvd_roc <- compute_roc(estimated_teststat_vec = -log10(esvd_pvalue),
                            true_de_idx = true_de_idx)
    esvd_roc <- smooth_roc(tpr = esvd_roc$tpr,
                           fpr = esvd_roc$fpr)
    esvd_auc <- compute_auc(tpr = esvd_roc$tpr,
                            fpr = esvd_roc$fpr)

    deseq_pvalue <- deseq2_res[paste0("gene-", 1:p), "pvalue"]
    deseq2_roc <- compute_roc(estimated_teststat_vec = -log10(deseq_pvalue),
                              true_de_idx = true_de_idx)
    deseq2_roc <- smooth_roc(tpr = deseq2_roc$tpr,
                             fpr = deseq2_roc$fpr)
    deseq2_auc <- compute_auc(tpr = deseq2_roc$tpr,
                              fpr = deseq2_roc$fpr)

    sctransform_pvalue <- de_result[paste0("gene-", 1:p), "p_val"]
    sctransform_roc <- compute_roc(estimated_teststat_vec = -log10(sctransform_pvalue),
                                   true_de_idx = true_de_idx)
    sctransform_roc <- smooth_roc(tpr = sctransform_roc$tpr,
                                  fpr = sctransform_roc$fpr)
    sctransform_auc <- compute_auc(tpr = sctransform_roc$tpr,
                                   fpr = sctransform_roc$fpr)

    mast_pvalue <- mast_pval_glmer[paste0("gene-", 1:p)]
    mast_roc <- compute_roc(estimated_teststat_vec = -log10(mast_pval_glmer),
                            true_de_idx = true_de_idx)
    mast_roc <- smooth_roc(tpr = mast_roc$tpr,
                           fpr = mast_roc$fpr)
    mast_auc <- compute_auc(tpr = mast_roc$tpr,
                            fpr = mast_roc$fpr)

    auc_vec <- c(esvd = esvd_auc,
                 deseq2 = deseq2_auc,
                 mast = mast_auc,
                 sct = sctransform_auc)

    print(round(auc_vec*100, 1))
  }
}



