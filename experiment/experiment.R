rm(list=ls())

load("../../../out/Writeup11e/Writeup11e_sns_invip_esvd9.RData")
eSVD_obj2 <- eSVD_obj
rm(list = "eSVD_obj")
load("../../../out/main/sns_invip_esvd.RData")

head(eSVD_obj$covariates[,c("Log_UMI", "age", "RNA.Integrity.Number", "post.mortem.hours", "percent.mt", "nFeature_RNA")])
head(eSVD_obj2$covariates[,c("Log_UMI", "age", "RNA.Integrity.Number", "post.mortem.hours", "percent.mt", "nFeature_RNA")])

col_names <- colnames(eSVD_obj$covariates)[which(colnames(eSVD_obj$covariates) %in% colnames(eSVD_obj2$covariates))]
sum(abs(eSVD_obj$covariates[,col_names] - eSVD_obj2$covariates[,col_names]))

for(col_name in col_names){
  print(col_name)
  print(sum(abs(eSVD_obj$covariates[,col_name] - eSVD_obj2$covariates[,col_name])))
  print("====")
}

######################

covariates2 <- eSVD_obj$covariates
