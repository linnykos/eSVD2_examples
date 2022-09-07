png("../../out/fig/simulation/tmp.png",
    height = 2500, width = 2500,
    units = "px", res = 300)
image(nat_mat)
graphics.off()

png("../../out/fig/simulation/tmp.png",
    height = 2500, width = 2500,
    units = "px", res = 300)
image(obs_mat)
graphics.off()

png("../../out/fig/simulation/tmp.png",
    height = 2500, width = 2500,
    units = "px", res = 300)
image(log1p(obs_mat))
graphics.off()
