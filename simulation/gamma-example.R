rm(list=ls())

set.seed(10)
n_subjects <- 10
case_color_palette <- grDevices::colorRampPalette(c(rgb(140, 0, 0, maxColorValue = 255),
                                                    rgb(244, 84, 84, maxColorValue = 255)))(n_subjects)
control_color_palette <- grDevices::colorRampPalette(c(rgb(47, 60, 190, maxColorValue = 255),
                                                       rgb(27, 198, 245, maxColorValue = 255)))(n_subjects)
transparent_gray <- rgb(0.5,0.5,0.5,0.4)
two_letters <- substr(transparent_gray, start = 8, stop = 9)
case_color_trans_palette <- paste0(case_color_palette, two_letters)
control_color_trans_palette <- paste0(control_color_palette, two_letters)

case_shape <- 36
case_rate <- 3
control_shape <- 3
control_rate <- 1/4

case_shape_vec <- case_shape + stats::runif(n_subjects, min = -1, max = 1)
case_rate_vec <- case_rate + stats::runif(n_subjects, min = -0.1, max = 0.1)
control_shape_vec <- control_shape + stats::runif(n_subjects, min = -0.05, max = 0.05)
control_rate_vec <- control_rate + stats::runif(n_subjects, min = -0.05, max = 0.05)

xlim <- c(0,40); x_adj <- 4
case_gaussian_list <- lapply(1:n_subjects, function(i){
  shape_val <- case_shape_vec[i]
  rate_val <- case_rate_vec[i]
  xseq <- seq(xlim[1], xlim[2], length.out = 1000)
  xseq <- pmax(xseq, 0)
  yseq <- stats::dgamma(xseq, shape = shape_val, rate = rate_val)
  yseq <- yseq - min(yseq)
  yseq <- yseq/max(yseq)*4
  cbind(xseq/x_adj, yseq)
})

control_gaussian_list <- lapply(1:n_subjects, function(i){
  shape_val <- control_shape_vec[i]
  rate_val <- control_rate_vec[i]
  xseq <- seq(xlim[1], xlim[2], length.out = 1000)
  xseq <- pmax(xseq, 0)
  yseq <- stats::dgamma(xseq, shape = shape_val, rate = rate_val)
  yseq <- yseq - min(yseq)
  yseq <- yseq/max(yseq)*4
  cbind(xseq/x_adj, yseq)
})


png(paste0("../../out/fig/simulation/gamma_within-indiv.png"),
    height = 800, width = 1500,
    units = "px", res = 500)
par(mar = c(3,0,0,0), bg = NA)
plot(NA,
     xlim = xlim/x_adj,
     ylim = c(0, 5.1),
     xaxt = "n", yaxt = "n", bty = "n",
     xlab = "", ylab = "")
for(i in 1:n_subjects){
  graphics::polygon(x = c(control_gaussian_list[[i]][,1], rev(control_gaussian_list[[i]][,1])),
                    y = 1+c(control_gaussian_list[[i]][,2], rep(0, nrow(control_gaussian_list[[i]]))),
                    col = control_color_trans_palette[i])
}
for(i in 1:n_subjects){
  graphics::polygon(x = c(case_gaussian_list[[i]][,1], rev(case_gaussian_list[[i]][,1])),
                    y = c(case_gaussian_list[[i]][,2], rep(0, nrow(case_gaussian_list[[i]]))),
                    col = case_color_trans_palette[i])
}

# for(i in 1:n_subjects){
#   lines(rep(control_shape_vec[i]/control_rate_vec[i]/x_adj, 2), c(0,10), lwd = 3, col = "white")
#   lines(rep(control_shape_vec[i]/control_rate_vec[i]/x_adj, 2), c(0,10), lwd = 2, lty = 2, col = control_color_palette[i])
#   lines(rep(case_shape_vec[i]/case_rate_vec[i]/x_adj, 2), c(0,10), lwd = 3, col = "white")
#   lines(rep(case_shape_vec[i]/case_rate_vec[i]/x_adj, 2), c(0,10), lwd = 2, lty = 2, col = case_color_palette[i])
# }

axis(side = 1, lwd = 2, cex.axis = 1.25, cex.lab = 1.25)
graphics.off()

##########################3

set.seed(10)
n_subjects <- 10
case_color_palette <- grDevices::colorRampPalette(c(rgb(140, 0, 0, maxColorValue = 255),
                                                    rgb(244, 84, 84, maxColorValue = 255)))(n_subjects)
control_color_palette <- grDevices::colorRampPalette(c(rgb(47, 60, 190, maxColorValue = 255),
                                                       rgb(27, 198, 245, maxColorValue = 255)))(n_subjects)
transparent_gray <- rgb(0.5,0.5,0.5,0.4)
two_letters <- substr(transparent_gray, start = 8, stop = 9)
case_color_trans_palette <- paste0(case_color_palette, two_letters)
control_color_trans_palette <- paste0(control_color_palette, two_letters)

case_shape <- 36
case_rate <- 3
control_shape_2 <- 25
control_rate_2 <- 1
control_rate_1 <- 1
control_shape_1 <- (((120-control_shape_2)/(n_subjects-1)))*control_rate_1

case_shape_vec <- case_shape + stats::runif(n_subjects, min = -15, max = 15)
case_rate_vec <- rep(case_rate, n_subjects)
control_shape_vec <- c(rep(control_shape_1, n_subjects-1), control_shape_2) + stats::runif(n_subjects, min = -5, max = 5)
control_rate_vec <- c(rep(control_rate_1, n_subjects-1), control_rate_2)

xlim <- c(0,40); x_adj <- 4
case_gaussian_list <- lapply(1:n_subjects, function(i){
  shape_val <- case_shape_vec[i]
  rate_val <- case_rate_vec[i]
  xseq <- seq(xlim[1], xlim[2], length.out = 1000)
  xseq <- pmax(xseq, 0)
  yseq <- stats::dgamma(xseq, shape = shape_val, rate = rate_val)
  yseq <- yseq - min(yseq)
  yseq <- yseq/max(yseq)*4
  cbind(xseq/x_adj, yseq)
})

control_gaussian_list <- lapply(1:n_subjects, function(i){
  shape_val <- control_shape_vec[i]
  rate_val <- control_rate_vec[i]
  xseq <- seq(xlim[1], xlim[2], length.out = 1000)
  xseq <- pmax(xseq, 0)
  yseq <- stats::dgamma(xseq, shape = shape_val, rate = rate_val)
  yseq <- yseq - min(yseq)
  yseq <- yseq/max(yseq)*4
  cbind(xseq/x_adj, yseq)
})


png(paste0("../../out/fig/simulation/gamma_between-indiv.png"),
    height = 800, width = 1500,
    units = "px", res = 500)
par(mar = c(3,0,0,0), bg = NA)
plot(NA,
     xlim = xlim/x_adj,
     ylim = c(0, 5.1),
     xaxt = "n", yaxt = "n", bty = "n",
     xlab = "", ylab = "")
for(i in 1:n_subjects){
  graphics::polygon(x = c(control_gaussian_list[[i]][,1], rev(control_gaussian_list[[i]][,1])),
                    y = 1+c(control_gaussian_list[[i]][,2], rep(0, nrow(control_gaussian_list[[i]]))),
                    col = control_color_trans_palette[i])
}
for(i in 1:n_subjects){
  graphics::polygon(x = c(case_gaussian_list[[i]][,1], rev(case_gaussian_list[[i]][,1])),
                    y = c(case_gaussian_list[[i]][,2], rep(0, nrow(case_gaussian_list[[i]]))),
                    col = case_color_trans_palette[i])
}

axis(side = 1, lwd = 2, cex.axis = 1.25, cex.lab = 1.25)
graphics.off()
