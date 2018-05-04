rm(list=ls())

#### DUMMY ####
load("../simulation/comparison_pvalue_knownsigma.RData")
zz <- res
load("../simulation/comparison_pvalue_unknownsigma.RData")
for(i in 1:length(res)){
  res[[i]][2,] <- res[[i]][1,]
  res[[i]][1,] <- zz[[i]][1,]
}

pdf("../figures/qq.pdf", height = 3, width = 10)
par(mfrow = c(1, length(res)))
for(i in 1:length(res)){
  plot(sort(res[[i]][1,]), seq(0, 1, length.out = length(res[[i]][1,])), xlim = c(0,1),
       ylim = c(0,1), pch = 16, asp = T, main = paste0("Jump size: ", paramMat[i,1]),
       xlab = "Observed quantile", ylab = "Theoretical quantile")
  points(sort(res[[i]][2,]), seq(0, 1, length.out = length(res[[i]][2,])), pch = 16, col = "red")
  legend("bottomright", c("Known sigma", "Unknown sigma"), cex=0.8,
         fill=c("black","red"), bty = "n")
  lines(c(0,1), c(0,1), col = "blue", lwd = 2)
}
graphics.off()

pdf("../figures/pairwise.pdf", height = 3, width = 10)
par(mfrow = c(1, length(res)))
for(i in 1:length(res)){
  plot(res[[i]][1,], res[[i]][2,], pch = 16, xlab = "Known sigma", ylab = "Unknown sigma",
       main = paste0("Jump size: ", paramMat[i,1]), col = rgb(0.5, 0.5, 0.5, 0.5 ), cex = 2)
  lines(c(0,1), c(0,1), col = "blue", lwd = 2)
}
graphics.off()

par(mfrow = c(1, length(res)))
for(i in 1:length(res)){
  plot(log(res[[i]][1,]+1e-6), log(res[[i]][2,]+1e-6), pch = 16, xlab = "Saturated", ylab = "Selected", main = names(res)[i], col = rgb(0.5, 0.5, 0.5, 0.5 ), cex = 2)
}

par(mfrow = c(1, length(res)))
for(i in 1:length(res)){
  plot(sort(res[[i]][1,]), seq(0, 1, length.out = length(res[[i]][1,])), xlim = c(0,1),
       ylim = c(0,1), pch = 16, asp = T)
  lines(c(0,1), c(0,1), col = "blue", lwd = 2)
}
