rm(list=ls())
library(Seurat)
library(SummarizedExperiment)
library(DESeq2)

for(ii in 1:3){
  for(jj in 1:3){
    print(paste0("Gene setting: ", ii, ", Individual setting: ", jj))

    load(paste0("~/kzlinlab/projects/eSVD2/out/simulation/simulation-power_geneSetting",
                ii, "_individualSetting", jj,
                ".RData"))

    set.seed(10)

    seurat_obj$cc <- factor(seurat_obj$cc)
    seurat_obj$gender <- factor(seurat_obj$gender)
    seurat_obj$tobacco <- factor(seurat_obj$tobacco)
    seurat_obj <- Seurat::SCTransform(seurat_obj, method = "glmGamPoi",
                                      residual.features = rownames(seurat_obj[["RNA"]]),
                                      vars.to.regress = c("age", "gender", "tobacco"),
                                      verbose = T)
    Seurat::Idents(seurat_obj) <- "cc"

    Seurat::DefaultAssay(seurat_obj) <- "SCT"
    de_result <- Seurat::FindMarkers(seurat_obj, ident.1 = "1", ident.2 = "0",
                                     slot = "scale.data",
                                     test.use = "wilcox",
                                     logfc.threshold = 0,
                                     min.pct = 0,
                                     verbose = T)

    date_of_run <- Sys.time()
    session_info <- devtools::session_info()
    save(de_result,
         date_of_run, session_info,
         file = paste0("~/kzlinlab/projects/eSVD2/out/simulation/simulation-power-sctransform_geneSetting",
                       ii, "_individualSetting", jj,
                       ".RData"))
  }
}

print("Done! :)")

