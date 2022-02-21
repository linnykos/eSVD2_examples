rm(list=ls())
load("../../../../out/Writeup10/Writeup10_sns_invip_processed2.RData")

Seurat::Idents(sns) <- "diagnosis"
set.seed(10)
de_res_nb <- Seurat::FindMarkers(sns,
                              features = sns[["RNA"]]@var.features,
                              ident.1 = "ASD",
                              ident.2 = "Control",
                              test.use = "negbinom",
                              logfc.threshold = 0,
                              min.pct = 0)
save(de_res_nb,
     file = "../../../../out/Writeup10/Writeup10_sns_invip_de.RData")
