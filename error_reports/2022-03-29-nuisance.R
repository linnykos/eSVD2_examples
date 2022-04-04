rm(list=ls())
load("../eSVD2_examples/error_reports/2022-03-29-nuisance.RData")

for(i in 1:length(nuisance_error_list$error)){
  print(i)
  testthat::expect_error(
    eSVD2:::gamma_rate(x = nuisance_error_list$error[[i]][,"x"],
                       mu = nuisance_error_list$error[[i]][,"mu"],
                       s = nuisance_error_list$error[[i]][,"s"])
  )
}

####

for(i in 1:length(nuisance_error_list$nan)){
  print(i)
  val <-  eSVD2:::gamma_rate(x = nuisance_error_list$nan[[i]][,"x"],
                             mu = nuisance_error_list$nan[[i]][,"mu"],
                             s = nuisance_error_list$nan[[i]][,"s"])
  stopifnot(is.nan(val))
}

quantile(nuisance_error_list$nan[[1]][,"x"])
quantile(nuisance_error_list$nan[[1]][,"mu"])
quantile(nuisance_error_list$nan[[1]][,"s"])

####

for(i in 1:length(nuisance_error_list$large)){
  print(i)
  val <-  eSVD2:::gamma_rate(x = nuisance_error_list$large[[i]][,"x"],
                             mu = nuisance_error_list$large[[i]][,"mu"],
                             s = nuisance_error_list$large[[i]][,"s"])
  stopifnot(val > 1000)
}

####

for(i in 1:length(nuisance_error_list$normal)){
  val <-  eSVD2:::gamma_rate(x = nuisance_error_list$normal[[i]][,"x"],
                             mu = nuisance_error_list$normal[[i]][,"mu"],
                             s = nuisance_error_list$normal[[i]][,"s"])
  print(round(val, 2))
}

