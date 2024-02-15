rm(list=ls())
library(Seurat)
library(SummarizedExperiment)
library(DESeq2)

for(ii in 1:3){
  print(paste0("Working on gene setting: ", ii))

  for(jj in 1:3){
    print(paste0("Working on size_factor setting: ", jj))

    load(paste0("~/kzlinlab/projects/eSVD2/out/simulation/simulation-power_geneSetting",
                ii, "_individualSetting", jj,
                ".RData"))

    #######

    p <- length(true_fdr_vec)
    true_de_idx <- which(true_fdr_vec < 0.05)

    lfc_vec <- .compute_true_lfc(
      case_individuals = case_individuals,
      control_individuals = control_individuals,
      covariates = covariates,
      individual_vec = individual_vec,
      x_mat = x_mat,
      y_mat = y_mat,
      z_mat = z_mat
    )

    break_vec <- seq(min(lfc_vec), max(lfc_vec), length.out = 25)

    png(paste0("~/kzlinlab/projects/eSVD2/out/fig/simulation_geneSetting",
               ii, "_individualSetting", jj, "_lfc.png"),
        height = 1200, width = 1800,
        units = "px", res = 500)
    par(mar = c(3,3,0.4,0.1))
    graphics::hist(
      lfc_vec,
      breaks = break_vec,
      col = "gray",
      xlab = "True Log-fold change",
      ylab = "Frequency",
      main = ""
    )

    graphics::hist(
      lfc_vec[true_de_idx],
      breaks = break_vec,
      col = "coral1",
      add = T
    )

    graphics.off()
  }
}
