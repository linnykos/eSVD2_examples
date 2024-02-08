rm(list=ls())
library(Seurat)
library(SummarizedExperiment)
library(DESeq2)
library(UniIsoRegression)

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

    #######

    p <- length(true_fdr_vec)
    true_de_idx <- which(true_fdr_vec < 0.05)

    eSVD_obj <- eSVD2:::compute_pvalue(input_obj = eSVD_obj)
    esvd_pvalue <- 10^(-eSVD_obj$pvalue_list$log10pvalue[paste0("gene-", 1:p)])
    esvd_roc <- compute_roc(estimated_teststat_vec = -log10(esvd_pvalue),
                            true_de_idx = true_de_idx)
    esvd_roc <- smooth_roc(tpr = esvd_roc$tpr,
                           fpr = esvd_roc$fpr)
    esvd_point <- roc_fdr_point(pvalue_vec = esvd_pvalue,
                                true_de_idx = true_de_idx,
                                fdr_threshold = 0.1)

    deseq_pvalue <- deseq2_res[paste0("gene-", 1:p), "pvalue"]
    deseq2_roc <- compute_roc(estimated_teststat_vec = -log10(deseq_pvalue),
                              true_de_idx = true_de_idx)
    deseq2_roc <- smooth_roc(tpr = deseq2_roc$tpr,
                             fpr = deseq2_roc$fpr)
    deseq2_point <- roc_fdr_point(pvalue_vec = deseq_pvalue,
                                  true_de_idx = true_de_idx,
                                  fdr_threshold = 0.1)

    sctransform_pvalue <- de_result[paste0("gene-", 1:p), "p_val"]
    sctransform_roc <- compute_roc(estimated_teststat_vec = -log10(sctransform_pvalue),
                                   true_de_idx = true_de_idx)
    sctransform_roc <- smooth_roc(tpr = sctransform_roc$tpr,
                                  fpr = sctransform_roc$fpr)
    sctransform_point <- roc_fdr_point(pvalue_vec = sctransform_pvalue,
                                       true_de_idx = true_de_idx,
                                       fdr_threshold = 0.1)


    orange_col <- rgb(235, 134, 47, maxColorValue = 255)
    yellow_col <- rgb(255, 205, 114, maxColorValue = 255)
    blue_col <- rgb(48, 174, 255, maxColorValue = 255)
    purple_col <- rgb(122, 49, 126, maxColorValue = 255)

    png(paste0("~/kzlinlab/projects/eSVD2/out/fig/simulation_geneSetting",
               ii, "_individualSetting", jj, "_roc.png"),
        height = 1750, width = 1250,
        units = "px", res = 500)
    par(mar = c(3,3,0.4,0.1))
    plot(NA, xlim = c(0,1), ylim = c(0,1), asp = T,
         xaxt = "n", yaxt = "n", bty = "n",
         cex.lab = 1.25, type = "n",
         xlab = "", ylab = "")
    lines(sctransform_roc$fpr, sctransform_roc$tpr, col = blue_col, lwd = 4)
    lines(deseq2_roc$fpr, deseq2_roc$tpr, col = yellow_col, lwd = 4)
    lines(esvd_roc$fpr, esvd_roc$tpr, col = orange_col, lwd = 4)

    lines(c(0,1), c(0,1), col = 1, lwd = 2, lty = 2)
    axis(1, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
    axis(2, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
    graphics.off()
  }
}



