rm(list=ls())
library(Seurat)
library(MAST)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

load("../../../../data/sns_autism/sns_formatted.RData")
head(sns@meta.data)
keep_vec <- rep(0, ncol(sns))
keep_vec[which(sns@meta.data$celltype == "L2/3")] <- 1
sns[["keep"]] <- keep_vec
sns <- subset(sns, keep == 1)

# following the analysis in https://github.com/himelmallick/BenchmarkSingleCell/blob/master/Library/run_MAST.R
# and https://www.bioconductor.org/packages/release/bioc/vignettes/MAST/inst/doc/MAITAnalysis.html
categorical_var <- c("diagnosis", "individual", "region", "sex", "Capbatch", "Seqbatch")
numerical_var <- c("age", "RNA.Integrity.Number", "post.mortem.hours", "percent.mt", "nFeature_RNA")
metadata <- sns@meta.data[,c(categorical_var, numerical_var)]
for(var in categorical_var){
  metadata[,var] <- as.factor( metadata[,var])
}
metadata[,"diagnosis"] <- stats::relevel(metadata[,"diagnosis"], "Control")

mat <- as.matrix(sns[["RNA"]]@counts)
tpms <- log1p(mat)

# create the SingleCellAssay (sca) object. See https://www.rdocumentation.org/packages/MAST/versions/0.931/topics/SingleCellAssay
sca <- MAST::FromMatrix(exprsArray = tpms, 
                        cData = data.frame(wellKey = rownames(metadata), grp = metadata), 
                        fData = data.frame(primerid = rownames(sns[["RNA"]]@counts)))

# put ngeneson into the sca object
ngeneson <- apply(mat, 2, function(x) mean(x>0))
CD <- colData(sca)
CD$ngeneson <- ngeneson
CD$cngeneson <- CD$ngeneson-mean(ngeneson)
colData(sca) <- CD

set.seed(10)
zlmdata <- MAST::zlm(~ cngeneson + grp.diagnosis + grp.individual + grp.region + grp.sex + grp.Capbatch + grp.Seqbatch + grp.age + grp.RNA.Integrity.Number + grp.post.mortem.hours + grp.percent.mt + grp.nFeature_RNA, sca = sca)

save(zlmdata, 
     file = "../../../../out/writeup8c/writeup8c_sns_layer23_de_mast_zlm.RData")

set.seed(10)
fit <- MAST::lrTest(zlmdata, "grp.diagnosisASD")

save(zlmdata, fit,
     file = "../../../../out/writeup8c/writeup8c_sns_layer23_de_mast_zlm.RData")

