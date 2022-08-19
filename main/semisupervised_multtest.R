## see alg 1 in https://arxiv.org/pdf/2106.13501v2.pdf

semisupervised_multtest <- function(teststat_vec,
                                    null_idx,
                                    alpha,
                                    verbose = F){
  stopifnot(length(null_idx) > 0, all(null_idx > 0), all(null_idx <= length(teststat_vec)),
            all(null_idx %% 1 == 0), length(null_idx) < length(teststat_vec),
            alpha >= 0, alpha <= 1,
            length(names(teststat_vec)) == length(teststat_vec))

  y_vec <- teststat_vec[null_idx]
  x_vec <- teststat_vec[-null_idx]
  n <- length(y_vec); m <- length(x_vec)
  order_vec <- order(abs(c(y_vec, x_vec)), decreasing = T)
  # (order_vec <= n)[1:100]
  orderx_vec <- order(abs(x_vec), decreasing = T)

  iter <- 0
  fdp <- 1; v <- n; l <- m+n; k <- m
  min_fdp <- 1; min_k <- NULL

  while(fdp > alpha & k >= 1){
    if(verbose) print(paste0("Iteration ", iter, ": FDP of ", fdp))
    l <- l-1
    if(l == 0) break()
    if(order_vec[l] <= n) v <- v-1 else k <- k-1

    if(k == 0) fdp <- 1 else fdp <- ((v+1)/k) * (m/(n+1))
    if(fdp <= min_fdp){min_fdp <- fdp; min_k <- k}
    iter <- iter + 1
  }
  if(fdp <= alpha) min_k <- k

  x_cutoff <- x_vec[orderx_vec[min_k]]
  names(x_vec)[abs(x_vec) >= x_cutoff]
}
