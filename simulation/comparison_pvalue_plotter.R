rm(list=ls())
load("../simulation/comparison_pvalue.RData")

for(i in 1:3){
  dat <- res[[i]]
  tmp <- dat[,which(!apply(dat,2,function(x){any(is.na(x))}))]

  plot(sort(tmp[1,]), seq(0, 1, length.out = length(tmp[1,])), pch = 16, asp = T, main = names(res)[i])
  points(sort(tmp[2,]), seq(0, 1, length.out = length(tmp[2,])), pch = 16, col = "red")
  lines(c(0,1), c(0,1), col = "blue", lwd = 2)
}
