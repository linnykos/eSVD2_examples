rm(list=ls())
library(Seurat)
library(eSVD2)

load("../../../out/main/regevEpi_ta1_preprocessed.RData")
# table(regevEpi$Subject, regevEpi$Sample_Health)
# table(regevEpi$Subject, regevEpi$Subject_Disease)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

keep_vec <- rep(1, ncol(regevEpi))
keep_vec[which(regevEpi$Sample_Health == "Inflamed")] <- 0
regevEpi$keep <- keep_vec
regevEpi <- subset(regevEpi, keep == 1)

# take only half of the healthy subjects
tab <- table(regevEpi$Subject, regevEpi$Subject_Disease)
healthy_subj <- rownames(tab[tab[,"HC"] != 0,])
set.seed(10)
split1 <- sample(healthy_subj, size = round(length(healthy_subj)/2), replace = F)
split2 <- setdiff(healthy_subj, split1)
keep_vec <- rep(1, ncol(regevEpi))
# Non-inflamed analysis uses split1, Inflamed uses split2
if(any(regevEpi$Sample_Health == "Non-inflamed")){
  keep_vec[which(regevEpi$Subject %in% split2)] <- 0
} else {
  keep_vec[which(regevEpi$Subject %in% split1)] <- 0
}
regevEpi$keep <- keep_vec
regevEpi <- subset(regevEpi, keep == 1)

regevEpi[["percent.mt"]] <- Seurat::PercentageFeatureSet(regevEpi, pattern = "^MT-")


#########################

regevEpi$Subject_Gender <- factor(regevEpi$Subject_Gender)
regevEpi$Subject_Smoking <- factor(regevEpi$Subject_Smoking)
regevEpi$Subject_Location <- factor(regevEpi$Subject_Location)
set.seed(10)
regevEpi <- Seurat::SCTransform(regevEpi, method = "glmGamPoi",
                                residual.features = regevEpi[["RNA"]]@var.features,
                                vars.to.regress = c("Subject_Gender", "percent.mt", "Subject_Location", "Subject_Smoking"),
                                verbose = T)
Seurat::Idents(regevEpi) <- "Subject_Disease"
levels(regevEpi)

Seurat::DefaultAssay(regevEpi) <- "SCT"
de_result <- Seurat::FindMarkers(regevEpi, ident.1 = "Colitis", ident.2 = "HC",
                                 slot = "scale.data",
                                 test.use = "wilcox",
                                 logfc.threshold = 0,
                                 min.pct = 0,
                                 verbose = T)

save(regevEpi, de_result,
     date_of_run, session_info,
     file = "../../../out/main/regev_epi_ta1-noninflamed_sctransform.RData")

