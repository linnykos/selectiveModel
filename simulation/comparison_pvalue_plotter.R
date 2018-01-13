rm(list=ls())
load("../simulation/comparison_pvalue.RData")

par(mfrow = c(1, length(res)))
for(i in 1:length(res)){
  dat <- res[[i]]
  tmp <- dat[,which(!apply(dat,2,function(x){any(is.na(x))}))]

  plot(sort(tmp[1,]), seq(0, 1, length.out = length(tmp[1,])), xlim = c(0,1), ylim = c(0,1), pch = 16, asp = T, main = names(res)[i])
  points(sort(tmp[2,]), seq(0, 1, length.out = length(tmp[2,])), pch = 16, col = "red")
  lines(c(0,1), c(0,1), col = "blue", lwd = 2)
}

par(mfrow = c(1, length(res)))
for(i in 1:length(res)){
  plot(res[[i]][1,], res[[i]][2,], pch = 16, xlab = "Saturated", ylab = "Selected", main = names(res)[i], col = rgb(0.5, 0.5, 0.5, 0.5 ), cex = 2)
}


par(mfrow = c(1, length(res)))
for(i in 1:length(res)){
  plot(log(res[[i]][1,]+1e-6), log(res[[i]][2,]+1e-6), pch = 16, xlab = "Saturated", ylab = "Selected", main = names(res)[i], col = rgb(0.5, 0.5, 0.5, 0.5 ), cex = 2)
}
