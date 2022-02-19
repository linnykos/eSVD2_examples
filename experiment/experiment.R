rm(list=ls())
load("../../../../out/Writeup10/Writeup10_sns_invip_processed2.RData")

mat <- as.matrix(Matrix::t(sns[["RNA"]]@counts[sns[["RNA"]]@var.features,]))

load("../../../../data/sns_autism/velmeshev_genes.RData")
tmp <- velmeshev_de_gene_df_list[[1]]
tmp <- tmp[which(tmp[,"Cell type"] == "IN-VIP"),]
de_gene_specific <- tmp[,"Gene name"]
de_gene_specific <- intersect(de_gene_specific, colnames(mat))

covariates <- sns@meta.data
covariates <- covariates[,c("nCount_RNA", "nFeature_RNA", "percent.mt",
                            "individual", "region", "age", "sex",
                            "RNA.Integrity.Number", "post.mortem.hours",
                            "diagnosis", "Seqbatch")]
df <- data.frame(covariates)
df[,"nCount_RNA"] <- Matrix::rowSums(mat)
df[,"individual"] <- as.factor(df[,"individual"])
df[,"region"] <- as.factor(df[,"region"])
df[,"diagnosis"] <- as.factor(df[,"diagnosis"])
df[,"sex"] <- as.factor(df[,"sex"])
df[,"Seqbatch"] <- as.factor(df[,"Seqbatch"])
df[,"nFeature_RNA"] <- scale(df[,"nFeature_RNA"], center = T, scale = T)
df[,"percent.mt"] <- scale(df[,"percent.mt"], center = T, scale = T)
df[,"RNA.Integrity.Number"] <- scale(df[,"RNA.Integrity.Number"], center = T, scale = T)
df[,"age"] <- scale(df[,"age"], center = T, scale = T)
df[,"post.mortem.hours"] <- scale(df[,"post.mortem.hours"], center = T, scale = T)

# var_idx <- which(colnames(mat) == "SAT2")
# var_idx <- which(colnames(mat) == "MT-CO2")
anova_list <- lapply(1:length(de_gene_specific), function(i){
  print(paste0(i, " of ", length(de_gene_specific)))

  var_idx <- which(colnames(mat) == de_gene_specific[i])
  df <- cbind(mat[,var_idx], df)
  colnames(df)[1] <- c("value")

  m1 <- lme4::glmer(value ~ (1|individual) + percent.mt + RNA.Integrity.Number + post.mortem.hours + (1|Seqbatch) + age + diagnosis + sex + region,
                    data = df,
                    family = stats::poisson(link = "log"),
                    offset = log(df[,"nCount_RNA"]))
  # summary(m1)
  m2 <- lme4::glmer(value ~ (1|individual) + percent.mt + RNA.Integrity.Number + post.mortem.hours + (1|Seqbatch) + age + sex + region,
                    data = df,
                    family = stats::poisson(link = "log"),
                    offset = log(df[,"nCount_RNA"]))

  stats::anova(m2, m1)
})

p_val <- sapply(anova_list, function(obj){
  obj["Pr(>Chisq)"][2,1]
})
floor(p_val*100)

############

mt_gene <- colnames(mat)[grep("^MT-", colnames(mat))]
anova_list2 <- lapply(1:length(mt_gene), function(i){
  print(paste0(i, " of ", length(mt_gene)))

  var_idx <- which(colnames(mat) == mt_gene[i])
  df <- cbind(mat[,var_idx], df)
  colnames(df)[1] <- c("value")

  m1 <- lme4::glmer(value ~ (1|individual) + percent.mt + RNA.Integrity.Number + post.mortem.hours + (1|Seqbatch) + age + diagnosis + sex + region,
                    data = df,
                    family = stats::poisson(link = "log"),
                    offset = log(df[,"nCount_RNA"]))
  # summary(m1)
  m2 <- lme4::glmer(value ~ (1|individual) + percent.mt + RNA.Integrity.Number + post.mortem.hours + (1|Seqbatch) + age + sex + region,
                    data = df,
                    family = stats::poisson(link = "log"),
                    offset = log(df[,"nCount_RNA"]))

  stats::anova(m2, m1)
})
p_val <- sapply(anova_list2, function(obj){
  obj["Pr(>Chisq)"][2,1]
})
floor(p_val*100)

