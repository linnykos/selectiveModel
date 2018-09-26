rm(list=ls())
load("../simulation/main_powercurve_knownsigma_tmp.RData")

sapply(res, dim)

# compute power
idx <- which(sapply(bs_res, length) > 1)
bs_vector <- lapply(bs_res[idx], function(x){
  x <- as.vector(x[3:4,])
  x <- x[!is.na(x)]
})

sapply(bs_vector, length)

bs_power_vector <- sapply(bs_vector, function(x){
  length(which(x <= 0.05/2))/length(x)
})
bs_power_vector
plot(paramMat_bs[,1], bs_power_vector, pch = 16, ylim = c(0,1),
     xlab = "Signal-to-noise ratio", ylab = "Conditional power over 2000 p-values",
     main = "Power for two-sided Bonferroni-corrected selected model test")
lines(paramMat_bs[,1], bs_power_vector)

######
# compute power
idx <- which(sapply(fl_res, length) > 1)
fl_vector <- lapply(fl_res[idx], function(x){
  x <- as.vector(x[3:4,])
  x <- x[!is.na(x)]
})

sapply(fl_vector, length)

fl_power_vector <- sapply(fl_vector, function(x){
  length(which(x <= 0.05/2))/length(x)
})
fl_power_vector
points(paramMat_fl[,1], fl_power_vector, pch = 16, col = "red")
lines(paramMat_fl[,1], fl_power_vector, col = "red")

legend("bottomright", c("2-step binary segementation",
                    "2-step fused lasso"), fill= c("black", "red"))


#########

.visualize_curve <- function(x, ...){
  pvec <- sapply(1:length(x), function(z){
    length(which(x[1:z] <= 0.05/2))/z
  })

  val <- pvec[length(pvec)]

  plot(pvec, ylim = c(0,1), ...)
  lines(c(-1e6, 1e6), rep(val, 2), col = "red", ...)

  invisible()
}

.visualize_curve(bs_vector[[1]])

.variance <- function(x){
  ind <- sapply(x, function(y){
    y <= 0.05/2
  })
  ind <- as.numeric(ind)
  stats::var(ind)/length(x)
}

.variance(bs_vector[[1]])
.variance(bs_vector[[1]][1:500])

####################################

i <- 1
vec <- as.numeric(res[[i]][9:12,])
vec <- vec[!is.na(vec)]
plot(sort(vec), seq(0, 1, length.out = length(vec)), asp = T)
lines(c(0,1), c(0,1), col = "red")
bool <- as.numeric(vec <= 0.05/2)

.bootstrap_power <- function(vec, trials = 1000){
  sd(sapply(1:trials, function(x){
    set.seed(x)
    vec2 <- sample(vec, length(vec), replace = T)
    mean(as.numeric(vec2 <= 0.05/2))
  }))
}

.bootstrap_power(vec)
