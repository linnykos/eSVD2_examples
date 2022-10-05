rm(list=ls())
case_color_palette <- grDevices::colorRampPalette(c(rgb(140, 0, 0, maxColorValue = 255),
                                                    rgb(244, 84, 84, maxColorValue = 255)))(3)
control_color_palette <- grDevices::colorRampPalette(c(rgb(47, 60, 190, maxColorValue = 255),
                                                       rgb(27, 198, 245, maxColorValue = 255)))(3)
transparent_gray <- rgb(0.5,0.5,0.5,0.4)
two_letters <- substr(transparent_gray, start = 8, stop = 9)
case_color_trans_palette <- paste0(case_color_palette, two_letters)
control_color_trans_palette <- paste0(control_color_palette, two_letters)

case_mean_vec <- c(0.9,1,1.1)
case_sd_vec <- c(0.8,1,0.9)
control_mean_vec <- c(1.5,1.6,1.7)
control_sd_vec <- c(0.8,1,0.9)

xlim <- c(0,4)
case_gaussian_list <- lapply(1:3, function(i){
  mean_val <- case_mean_vec[i]
  sd_val <- case_sd_vec[i]
  xseq <- seq(xlim[1], xlim[2], length.out = 1000)
  xseq <- pmax(xseq, 0)
  yseq <- stats::dnorm(xseq, mean = mean_val, sd = sd_val)
  yseq <- yseq - min(yseq)
  yseq <- yseq/max(yseq)*4
  cbind(xseq, yseq)
})
control_gaussian_list <- lapply(1:3, function(i){
  mean_val <- control_mean_vec[i]
  sd_val <- control_sd_vec[i]
  xseq <- seq(xlim[1], xlim[2], length.out = 1000)
  xseq <- pmax(xseq, 0)
  yseq <- stats::dnorm(xseq, mean = mean_val, sd = sd_val)
  yseq <- yseq - min(yseq)
  yseq <- yseq/max(yseq)*4
  cbind(xseq, yseq)
})


png(paste0("../../out/fig/simulation/gaussian_within-indiv.png"),
    height = 800, width = 1500,
    units = "px", res = 500)
par(mar = c(3,0,0,0), bg = NA)
plot(NA,
     xlim = xlim,
     ylim = c(0, 5.1),
     xaxt = "n", yaxt = "n", bty = "n",
     xlab = "", ylab = "")
for(i in 1:3){
  graphics::polygon(x = c(control_gaussian_list[[i]][,1], rev(control_gaussian_list[[i]][,1])),
                    y = 1+c(control_gaussian_list[[i]][,2], rep(0, nrow(control_gaussian_list[[i]]))),
                    col = control_color_trans_palette[i])
}
for(i in 1:3){
  graphics::polygon(x = c(case_gaussian_list[[i]][,1], rev(case_gaussian_list[[i]][,1])),
                    y = c(case_gaussian_list[[i]][,2], rep(0, nrow(case_gaussian_list[[i]]))),
                    col = case_color_trans_palette[i])
}

for(i in 1:3){
  lines(rep(control_mean_vec[i], 2), c(0,10), lwd = 3, col = "white")
  lines(rep(control_mean_vec[i], 2), c(0,10), lwd = 2, lty = 2, col = control_color_palette[i])
  lines(rep(case_mean_vec[i], 2), c(0,10), lwd = 3, col = "white")
  lines(rep(case_mean_vec[i], 2), c(0,10), lwd = 2, lty = 2, col = case_color_palette[i])
}

axis(side = 1, lwd = 2, cex.axis = 1.25, cex.lab = 1.25)
graphics.off()

##########################3


case_mean_vec <- c(0.5,1.5,3.8)
case_sd_vec <- c(0.2,0.5,0.2)
case_size_vec <-  c(0.75, 2, 0.75)
control_mean_vec <- c(2.1,4,4.5)
control_sd_vec <- c(0.2,0.5,0.2)
control_size_vec <- c(0.75, 2, 0.75)


xlim <- c(0,5.5)
case_gaussian_list <- lapply(1:3, function(i){
  mean_val <- case_mean_vec[i]
  sd_val <- case_sd_vec[i]
  xseq <- seq(xlim[1], xlim[2], length.out = 1000)
  xseq <- pmax(xseq, 0)
  yseq <- stats::dnorm(xseq, mean = mean_val, sd = sd_val)
  yseq <- yseq - min(yseq)
  yseq <- yseq/max(yseq)*2*case_size_vec[i]
  cbind(xseq, yseq)
})
control_gaussian_list <- lapply(1:3, function(i){
  mean_val <- control_mean_vec[i]
  sd_val <- control_sd_vec[i]
  xseq <- seq(xlim[1], xlim[2], length.out = 1000)
  xseq <- pmax(xseq, 0)
  yseq <- stats::dnorm(xseq, mean = mean_val, sd = sd_val)
  yseq <- yseq - min(yseq)
  yseq <- yseq/max(yseq)*2*control_size_vec[i]
  cbind(xseq, yseq)
})


png(paste0("../../out/fig/simulation/gaussian_between-indiv.png"),
    height = 800, width = 1500,
    units = "px", res = 500)
par(mar = c(3,0,0,0), bg = NA)
plot(NA,
     xlim = xlim,
     ylim = c(0, 5.1),
     xaxt = "n", yaxt = "n", bty = "n",
     xlab = "", ylab = "")

graphics::polygon(x = c(control_gaussian_list[[2]][,1], rev(control_gaussian_list[[2]][,1])),
                  y = 1+c(control_gaussian_list[[2]][,2], rep(0, nrow(control_gaussian_list[[2]]))),
                  col = control_color_palette[2])
for(i in c(1,3)){
  graphics::polygon(x = c(control_gaussian_list[[i]][,1], rev(control_gaussian_list[[i]][,1])),
                    y = 1+c(control_gaussian_list[[i]][,2], rep(0, nrow(control_gaussian_list[[i]]))),
                    col = control_color_palette[i])
}
graphics::polygon(x = c(case_gaussian_list[[2]][,1], rev(case_gaussian_list[[2]][,1])),
                  y = c(case_gaussian_list[[2]][,2], rep(0, nrow(case_gaussian_list[[2]]))),
                  col = case_color_palette[2])
for(i in c(1,3)){
  graphics::polygon(x = c(case_gaussian_list[[i]][,1], rev(case_gaussian_list[[i]][,1])),
                    y = c(case_gaussian_list[[i]][,2], rep(0, nrow(case_gaussian_list[[i]]))),
                    col = case_color_palette[i])
}

for(i in 1:3){
  lines(rep(control_mean_vec[i], 2), c(0,10), lwd = 3, col = "white")
  lines(rep(control_mean_vec[i], 2), c(0,10), lwd = 2, lty = 2, col = control_color_palette[i])
  lines(rep(case_mean_vec[i], 2), c(0,10), lwd = 3, col = "white")
  lines(rep(case_mean_vec[i], 2), c(0,10), lwd = 2, lty = 2, col = case_color_palette[i])
}

axis(side = 1, lwd = 2, cex.axis = 1.25, cex.lab = 1.25)
graphics.off()
