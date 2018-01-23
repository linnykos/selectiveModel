rm(list=ls())
load("../simulation/pvalue_variance.RData")
paramMat <- as.matrix(expand.grid(1:4, c(100, 1000, 7500), c(1, 10, 100), c(100, 500, 1000)))
colnames(paramMat) <- c("Level", "Burn.in", "Lapse", "Samples")

names(res) <- apply(paramMat, 1, paste0, collapse = "-")
mean_vec <- as.vector(sapply(res, mean))
var_vec <- as.vector(sapply(res, var))

#####

len <- length(unique(paramMat[,"Level"]))
mat <- vector("list", len)
for(i in 1:len){
  idx <- which(paramMat[,1] == i)
  mat[[i]] <- data.frame(name = names(res[idx]), mean = mean_vec[idx],
                         sd = sqrt(var_vec[idx]), lower = mean_vec[idx] - sqrt(var_vec[idx]),
                         upper = mean_vec[idx] + sqrt(var_vec[idx]))
}

#####

len <- length(unique(paramMat[,"Level"]))
vec1 <- numeric(0)
vec2 <- numeric(0)
for(i in 1:length(res)){
  vec1 <- c(vec1, res[[i]])
  vec2 <- c(vec2, rep(names(res)[i], length(res[[i]])))
}
dat <- data.frame(res = vec1, name = vec2)
ylim <- vector("list", 4)
for(i in 1:len){
  idx <- grep(paste0("^", i, "-"), vec2)
  ylim[[i]] <- c(min(vec1[idx]), max(vec1[idx]))
}

#####

pdf("../figures/pvalue_variance_base.pdf", 10, 8)
level_vec <- seq(0, 1.5, length.out = 4)
len <- length(unique(paramMat[,"Level"]))
par(mfrow = c(2,2))
for(i in 1:len){
  labels <- paste0(i, c("-100-1-100", "-1000-10-500", "-7500-100-1000"))
  idx <- which(dat$name %in% labels)
  dat2 <- dat[idx,]
  dat2[,2] <- factor(dat2[,2], levels = labels)
  boxplot(res ~ name, dat2, col = "gray", main = paste0("Level: ", level_vec[i]), pch = 16,
          ylim = ylim[[i]])
}
graphics.off()

######

pdf("../figures/pvalue_variance_burnin_small.pdf", 10.5, 8)
level_vec <- seq(0, 1.5, length.out = 4)
len <- length(unique(paramMat[,"Level"]))
par(mfrow = c(2,2))
for(i in 1:len){
  labels <- paste0(i, c("-100-1-100", "-1000-1-100", "-7500-1-100"))
  idx <- which(dat$name %in% labels)
  dat2 <- dat[idx,]
  dat2[,2] <- factor(dat2[,2], levels = labels)
  boxplot(res ~ name, dat2, col = "gray", main = paste0("Level: ", level_vec[i]), pch = 16,
          ylim = ylim[[i]])
}
graphics.off()

######

pdf("../figures/pvalue_variance_burnin_large.pdf", 10.5, 8)
level_vec <- seq(0, 1.5, length.out = 4)
len <- length(unique(paramMat[,"Level"]))
par(mfrow = c(2,2))
for(i in 1:len){
  labels <- paste0(i, c("-100-100-1000", "-1000-100-1000", "-7500-100-1000"))
  idx <- which(dat$name %in% labels)
  dat2 <- dat[idx,]
  dat2[,2] <- factor(dat2[,2], levels = labels)
  boxplot(res ~ name, dat2, col = "gray", main = paste0("Level: ", level_vec[i]), pch = 16,
          ylim = ylim[[i]])
}
graphics.off()

######

pdf("../figures/pvalue_variance_samples.pdf", 10.5, 8)
level_vec <- seq(0, 1.5, length.out = 4)
len <- length(unique(paramMat[,"Level"]))
par(mfrow = c(2,2))
for(i in 1:len){
  labels <- paste0(i, c("-100-1-100", "-100-1-500", "-100-1-1000"))
  idx <- which(dat$name %in% labels)
  dat2 <- dat[idx,]
  dat2[,2] <- factor(dat2[,2], levels = labels)
  boxplot(res ~ name, dat2, col = "gray", main = paste0("Level: ", level_vec[i]), pch = 16,
          ylim = ylim[[i]])
}
graphics.off()

######

pdf("../figures/pvalue_variance_lapse.pdf", 10.5, 8)
level_vec <- seq(0, 1.5, length.out = 4)
len <- length(unique(paramMat[,"Level"]))
par(mfrow = c(2,2))
for(i in 1:len){
  labels <- paste0(i, c("-100-1-100", "-100-10-100", "-100-100-100"))
  idx <- which(dat$name %in% labels)
  dat2 <- dat[idx,]
  dat2[,2] <- factor(dat2[,2], levels = labels)
  boxplot(res ~ name, dat2, col = "gray", main = paste0("Level: ", level_vec[i]), pch = 16,
          ylim = ylim[[i]])
}
graphics.off()

