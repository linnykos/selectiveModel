rm(list=ls())
#load("../simulation/comparison_pvalue.RData")
#load("../simulation/comparison_pvalue_signal2.RData")
load("../simulation/tmp.RData")

res <- res[[1]]
# vec1 <- res[1,]; vec1 <- vec1[!is.nan(vec1)]
# vec2 <- res[2,]; vec2 <- vec2[!is.nan(vec2)]
# min_val <- min(c(vec1, vec2), na.rm = T)
# max_val <- max(c(vec1, vec2), na.rm = T)

# par(mfrow = c(1,2))
# hist(vec1, xlim = c(min_val, max_val), breaks = 15, col = rgb(1,0,0,0.5))
# hist(vec2, breaks = 15, add = T, col = rgb(0,0,1,0.5))
#
tmp <- res[,which(!apply(res,2,function(x){any(is.na(x))}))]
# plot(tmp[1,], tmp[2,], xlab = "Saturated p value", ylab = "Selected p value", pch = 16)

par(mfrow = c(1,3))

plot(sort(tmp[1,]), seq(0, 1, length.out = length(tmp[1,])), pch = 16, asp = T)
points(sort(tmp[2,]), seq(0, 1, length.out = length(tmp[2,])), pch = 16, col = "red")
lines(c(0,1), c(0,1), col = "blue", lwd = 2)

####

load("../simulation/comparison_pvalue.RData")
res <- res[[1]]
tmp <- res[,which(!apply(res,2,function(x){any(is.na(x))}))]
plot(sort(tmp[1,]), seq(0, 1, length.out = length(tmp[1,])), pch = 16, asp = T)
points(sort(tmp[2,]), seq(0, 1, length.out = length(tmp[2,])), pch = 16, col = "red")
lines(c(0,1), c(0,1), col = "blue", lwd = 2)

load("../simulation/comparison_pvalue_signal2.RData")
res <- res[[1]]
tmp <- res[,which(!apply(res,2,function(x){any(is.na(x))}))]
plot(sort(tmp[1,]), seq(0, 1, length.out = length(tmp[1,])), pch = 16, asp = T)
points(sort(tmp[2,]), seq(0, 1, length.out = length(tmp[2,])), pch = 16, col = "red")
lines(c(0,1), c(0,1), col = "blue", lwd = 2)
