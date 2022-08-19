rm(list=ls())
library(Seurat)
library(eSVD2)

file_vec <- sort(c("regevEpi_ta1-inflamed_esvd.RData", "regevEpi_ta1-noninflamed_esvd.RData",
                   "regevEpi_ta2-inflamed_esvd.RData", "regevEpi_ta2-noninflamed_esvd.RData",
                   "regevEpi_entprog-inflamed_esvd.RData", "regevEpi_entprog-noninflamed_esvd.RData",
                   "regevEpi_cyclingta-inflamed_esvd.RData", "regevEpi_cyclingta-noninflamed_esvd.RData"))

tab_regevEpi_list <- sapply(1:length(file_vec), function(i){
  print(i)
  load(paste0("../../../out/main/", file_vec[i]))
  tab <- table(regevEpi$Sample, regevEpi$Subject_Disease)
  tab
})
sapply(tab_regevEpi_list, function(mat){
  colSums(mat)
})
sapply(tab_regevEpi_list, function(mat){
  case_vec <- mat[which(mat[,1]!=0),1]
  control_vec <-mat[which(mat[,2]!=0),2]

  print(paste0(length(case_vec), " - ", length(control_vec)))
  invisible()
})

file_vec <- sort(c("regevImm_macro-inflamed_esvd.RData", "regevImm_macro-noninflamed_esvd.RData",
                   "regevImm_plasma-inflamed_esvd.RData", "regevImm_plasma-noninflamed_esvd.RData"))
tab_regevImm_list <- sapply(1:length(file_vec), function(i){
  print(i)
  load(paste0("../../../out/main/", file_vec[i]))
  tab <- table(regevImm$Sample, regevImm$Subject_Disease)
  tab
})
sapply(tab_regevImm_list, function(mat){
  colSums(mat)
})
sapply(tab_regevImm_list, function(mat){
  case_vec <- mat[which(mat[,1]!=0),1]
  control_vec <-mat[which(mat[,2]!=0),2]

  print(paste0(length(case_vec), " - ", length(control_vec)))
  invisible()
})


file_vec <- c("adams_T_esvd.RData", "habermann_T_esvd.RData")
tab_lung_list <- sapply(1:length(file_vec), function(i){
  print(i)
  load(paste0("../../../out/main/", file_vec[i]))
  if(i == 1){
    tab <- table(adams$Subject_Identity, adams$Disease_Identity)
  } else {
    tab <- table(habermann$Sample_Name, habermann$Diagnosis)
  }
  tab
})
for(i in 1:length(tab_lung_list)){
  tab_lung_list[[i]] <- tab_lung_list[[i]][,c("IPF", "Control")]
}
sapply(tab_lung_list, function(mat){
  colSums(mat)
})
sapply(tab_lung_list, function(mat){
  case_vec <- mat[which(mat[,1]!=0),1]
  control_vec <-mat[which(mat[,2]!=0),2]

  print(paste0(length(case_vec), " - ", length(control_vec)))
  invisible()
})

#######################################

histogram_function <- function(tab_list, width, minor_gap, major_gap){
  len <- length(tab_list)
  par(mar = c(2,0,0.1,0))
  y_max <- width*2*len + minor_gap*(len) + major_gap*(len-1)
  plot(NA, xlim = c(0, log10(35000)),
       ylim = c(0, y_max),
       xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
  x_break_vals <- c(0,unlist(lapply(1:4, function(x){
    log10(seq(10^(x-1), 10^x, length.out = 6))[-1]
  })))
  for(x in x_break_vals){
    if(abs(x %% 1) <= 1e-4){
      lines(rep(x,2), c(-1e4,1e4), lwd = 2, lty = 1, col = "gray")
    } else {
      lines(rep(x,2), c(-1e4,1e4), lwd = 1, lty = 2, col = "gray")
    }
  }
  for(i in 1:len){
    case_vec <- tab_list[[i]][which(tab_list[[i]][,1]!=0),1]
    control_vec <- tab_list[[i]][which(tab_list[[i]][,2]!=0),2]

    case_total_len <- log10(sum(case_vec))
    control_total_len <- log10(sum(control_vec))

    num_case_subj <- length(case_vec)
    num_control_subj <- length(control_vec)

    base_palette <- RColorBrewer::brewer.pal(11, name = "RdYlBu")
    case_color_palette <- grDevices::colorRampPalette(base_palette[1:4])(num_case_subj)
    control_color_palette <- grDevices::colorRampPalette(base_palette[8:11])(num_control_subj)

    y_base <- (i-1)*(2*width + minor_gap + major_gap)
    lines(x = c(case_total_len, log10(35000)),
          y = rep(mean(c(y_base, y_base+width)), 2),
          col = 2, lwd = 2, lty = 1)
    lines(x = c(control_total_len, log10(35000)),
          y = rep(mean(c(y_base+width+minor_gap, y_base+2*width+minor_gap)), 2),
          col = 2, lwd = 2, lty = 1)

    polygon(x = c(0, case_total_len, case_total_len, 0),
            y = c(rep(y_base, 2), rep(y_base+width, 2)),
            col = "black", lwd = 3, border = "black")

    polygon(x = c(0, control_total_len, control_total_len, 0),
            y = c(rep(y_base+width+minor_gap, 2), rep(y_base+2*width+minor_gap, 2)),
            col = "black", lwd = 3, border = "black")

    case_prop_len <- sort((case_vec/sum(case_vec))*case_total_len, decreasing = T)
    case_prop_len_cum <- c(0,cumsum(case_prop_len))
    for(i in 1:num_case_subj){
      polygon(x = c(case_prop_len_cum[i], case_prop_len_cum[i+1], case_prop_len_cum[i+1], case_prop_len_cum[i]),
              y = c(rep(y_base, 2), rep(y_base+width, 2)),
              col = case_color_palette[i], border = "black")
    }

    control_prop_len <- sort((control_vec/sum(control_vec))*control_total_len, decreasing = T)
    control_prop_len_cum <- c(0,cumsum(control_prop_len))
    for(i in 1:num_control_subj){
      polygon(x = c(control_prop_len_cum[i], control_prop_len_cum[i+1], control_prop_len_cum[i+1], control_prop_len_cum[i]),
              y =c(rep(y_base+width+minor_gap, 2), rep(y_base+2*width+minor_gap, 2)),
              col = control_color_palette[i], border = "black")
    }

    axis(1, at = seq(0, 4, by = 1),
         labels = c("0", "10", "100", "1000", "10000"),
         cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
  }
}

width <- 1
minor_gap <- 0.5
major_gap <- 1.5

png("../../../out/fig/main/regevEpi_sample_histogram.png",
    height = 2700, width = 1500,
    units = "px", res = 500)
par(bg = NA)
histogram_function(tab_regevEpi_list,
                   width = width, minor_gap = minor_gap, major_gap = major_gap)
graphics.off()

png("../../../out/fig/main/regevImm_sample_histogram.png",
    height = 1350, width = 1500,
    units = "px", res = 500)
par(bg = NA)
histogram_function(tab_regevImm_list,
                   width = width, minor_gap = minor_gap, major_gap = major_gap)
graphics.off()

png("../../../out/fig/main/lung_sample_histogram.png",
    height = 720, width = 1500,
    units = "px", res = 500)
par(bg = NA)
histogram_function(tab_lung_list,
                   width = width, minor_gap = minor_gap, major_gap = major_gap)
graphics.off()

