rm(list=ls())
trials <- 1000
n1 <- 1000; n2 <- 1000
p_val <- sapply(1:trials, function(trial){
  set.seed(trial)
  x <- rnorm(n1)
  y <- rnorm(n2)

  mean1 <- mean(x); mean2 <- mean(y)
  cov1 <- var(x); cov2 <- var(y)
  combined_cov <- ((n1-1)*cov1 + (n2-1)*cov2)/(n1+n2-2)
  test_stat <- (n1*n2)*(mean1 - mean2)^2/(combined_cov * (n1+n2))

  test_stat

  p <- 1; m <- n1+n2-2
  p_val <- 1-stats::pf((m-p+1)/(p*m)*test_stat, df1 = p, df2 = m-p+1)

  p_val
})

plot(sort(p_val), seq(0,1,length.out=length(p_val)), asp = T)
