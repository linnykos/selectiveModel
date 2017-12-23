load("../simulation/four_jump.RData")

pdf("../figures/four_jump_selective_model.pdf", height = 5, width = 5)
plot(NA, xlim = c(0,1), ylim = c(0,1), main = "QQ-plot", ylab = "", xlab = "", asp = T)

for(i in 1:length(res)){
  vec <- res[[i]][!is.na(res[[i]])]
  points(sort(vec), seq(0,1,length.out = length(vec)), pch = 16, col = i)
}

lines(x = c(0,1), y = c(0,1), lwd = 5, lty = 2, col = "red")

legend("bottomright", c(paste0("Delta: ", paramMat[,1])),
       bty="n", fill=c(1:5));

graphics.off()

