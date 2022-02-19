initialize_lme4 <- function(dat,
                            k,
                            metadata,
                            pval_thres = 0.05,
                            tol = 1e-3,
                            verbose = 0,
                            tmp_path = NULL){
  df <- metadata[,c("percent.mt", "individual", "region", "age", "sex",
                              "RNA.Integrity.Number", "post.mortem.hours",
                              "diagnosis", "Seqbatch")]
  df <- data.frame(df)
  df[,"Log_UMI"] <- log(Matrix::rowSums(dat))
  df[,"individual"] <- as.factor(df[,"individual"])
  df[,"region"] <- as.factor(df[,"region"])
  df[,"diagnosis"] <- factor(df[,"diagnosis"], levels = c("Control", "ASD"))
  df[,"sex"] <- as.factor(df[,"sex"])
  df[,"Seqbatch"] <- as.factor(df[,"Seqbatch"])
  df[,"percent.mt"] <- scale(df[,"percent.mt"], center = T, scale = T)
  df[,"RNA.Integrity.Number"] <- scale(df[,"RNA.Integrity.Number"], center = T, scale = T)
  df[,"age"] <- scale(df[,"age"], center = T, scale = T)
  df[,"post.mortem.hours"] <- scale(df[,"post.mortem.hours"], center = T, scale = T)
  covariates <- .create_covariate_matrix(df)

  n <- nrow(dat); p <- ncol(dat)

  coef_mat <- matrix(0, nrow = p, ncol = ncol(covariates))
  rownames(coef_mat) <- colnames(dat)
  colnames(coef_mat) <- colnames(covariates)

  if(verbose >= 1) print("Starting step 1: Fitting lme4's to covariates")
  for(j in 1:p){
    if(verbose == 1 && p > 10 && j %% floor(p/10) == 0) cat('*')
    df_tmp <- cbind(dat[,j], df)
    colnames(df_tmp)[1] <- c("value")

    m1 <- lme4::glmer(value ~ (1|individual) + percent.mt + RNA.Integrity.Number + post.mortem.hours + (1|Seqbatch) + age + diagnosis + sex + region,
                      data = df_tmp,
                      family = stats::poisson(link = "log"),
                      offset = df_tmp[,"Log_UMI"],
                      control = lme4::glmerControl(check.conv.singular = lme4::.makeCC(action = "ignore",  tol = 1e-4)))

    m2 <- lme4::glmer(value ~ (1|individual) + percent.mt + RNA.Integrity.Number + post.mortem.hours + (1|Seqbatch) + age + sex + region,
                      data = df_tmp,
                      family = stats::poisson(link = "log"),
                      offset = df_tmp[,"Log_UMI"],
                      control = lme4::glmerControl(check.conv.singular = lme4::.makeCC(action = "ignore",  tol = 1e-4)))

    anova_res <- stats::anova(m2, m1)
    p_val <- anova_res["Pr(>Chisq)"]["m1","Pr(>Chisq)"]
    if(p_val < pval_thres){
      coef_vec <- .lme4_extract_coefficient(lme4_model = m1,
                                            covariate_names = colnames(covariates))
    } else {
      coef_vec <- .lme4_extract_coefficient(lme4_model = m2,
                                            covariate_names = colnames(covariates))
    }

    if(verbose == 2) print(paste0("Finished gene ", j, " (", colnames(dat)[j], "): Coefficient of ", coef_vec["diagnosis_ASD"]))
    coef_mat[j,] <- coef_vec
    if(!is.null(tmp_path) && p > 10 && j %% floor(p/10) == 0) {
      save(covariates, coef_mat, file = tmp_path)
    }
  }

  if(verbose >= 1) print("Starting step 2: Fitting SVD of residuals")
  dat_transform <- log1p(dat)
  nat_mat <- tcrossprod(covariates, coef_mat)
  residual_mat <- dat_transform - nat_mat

  svd_res <- eSVD2:::.svd_truncated(residual_mat,
                                    K = k,
                                    symmetric = F,
                                    rescale = F,
                                    mean_vec = NULL,
                                    sd_vec = NULL,
                                    K_full_rank = F)
  x_init <- eSVD2:::.mult_mat_vec(svd_res$u, sqrt(svd_res$d))
  y_init <- eSVD2:::.mult_mat_vec(svd_res$v, sqrt(svd_res$d))

  rownames(x_init) <- rownames(dat)
  rownames(y_init) <- colnames(dat)

  structure(list(x_mat = x_init, y_mat = y_init,
                 b_mat = coef_mat,
                 covariates = covariates,
                 nuisance_param_vec = rep(1, ncol(dat)),
                 offset_vec = rep(0, nrow(dat)),
                 pval_thres = pval_thres),
            class = "eSVD")
}

.create_covariate_matrix <- function(df){
  numerical_var <- c("Log_UMI", "percent.mt", "age", "RNA.Integrity.Number", "post.mortem.hours")
  twocat_category_var <- c("sex", "diagnosis", "region")
  multi_category_var <- c("individual", "Seqbatch")
  n <- nrow(df)

  covariates <- cbind(rep(1, n), df[,numerical_var])
  colnames(covariates) <- c("Intercept", numerical_var)
  for(variable in twocat_category_var){
    vec <- df[,variable]
    uniq_level <- levels(vec)[2]
    tmp <- rep(0, n)
    tmp[which(vec == uniq_level)] <- 1

    var_name <- paste0(variable, "_", uniq_level)
    covariates <- cbind(covariates, tmp)
    colnames(covariates)[ncol(covariates)] <- var_name
  }

  for(variable in multi_category_var){
    vec <- df[,variable]
    uniq_level <- sort(levels(vec)) #[[note to self: this line is hard-coded]]
    for(lvl in uniq_level){
      tmp <- rep(0, n)
      tmp[which(vec == lvl)] <- 1

      var_name <- paste0(variable, "_", lvl)
      covariates <- cbind(covariates, tmp)
      colnames(covariates)[ncol(covariates)] <- var_name
    }
  }

  as.matrix(covariates)
}

.lme4_extract_coefficient <- function(lme4_model,
                                      covariate_names){
  coef_vec <- rep(0, length(covariate_names))
  names(coef_vec) <- covariate_names

  lme4_vec1 <- stats::coef(summary(lme4_model))[,"Estimate"]
  coef_vec["Intercept"] <- lme4_vec1["(Intercept)"]
  coef_vec["Log_UMI"] <- 1
  numerical_var <- c("percent.mt", "age", "RNA.Integrity.Number", "post.mortem.hours")
  for(var in numerical_var){
    coef_vec[var] <- lme4_vec1[var]
  }

  twocat_category_var <- c("sex", "diagnosis", "region")
  for(var in twocat_category_var){
    idx1 <- grep(paste0("^", var), names(coef_vec))
    idx2 <- grep(paste0("^", var), names(lme4_vec1))
    if(length(idx2) > 0){
      coef_vec[idx1] <- lme4_vec1[idx2]
    } else {
      coef_vec[idx1] <- 0
    }
  }

  lme4_vec2 <- lme4::ranef(lme4_model)[["individual"]]
  # cbind(lme4_vec2+coef_vec[1], coef(lme4_model)$individual[,"(Intercept)"])
  for(indiv in rownames(lme4_vec2)){
    idx <- which(covariate_names == paste0("individual_", indiv))
    coef_vec[idx] <- lme4_vec2[indiv,1]
  }

  lme4_vec3 <- lme4::ranef(lme4_model)[["Seqbatch"]]
  for(batch in rownames(lme4_vec3)){
    idx <- which(covariate_names == paste0("Seqbatch_", batch))
    coef_vec[idx] <- lme4_vec3[batch,1]
  }

  coef_vec
}
