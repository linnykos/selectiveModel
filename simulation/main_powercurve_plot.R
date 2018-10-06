rm(list=ls())
load("../simulation/main_powercurve_unknownsigma.RData")
bs_res_unknown <- bs_res
fl_res_unknown <- fl_res

load("../simulation/main_powercurve_knownsigma.RData")
bs_res_known <- bs_res
fl_res_known <- fl_res

# compute power
.bootstrap_power <- function(vec, trials = 1000){
  sd(sapply(1:trials, function(x){
    set.seed(x)
    vec2 <- sample(vec, length(vec), replace = T)
    mean(as.numeric(vec2 <= 0.05/2))
  }))
}

.compute_power <- function(res1){
  if(is.list(res1)){
    idx <- sapply(res1, length)
    res1 <- res1[-which(idx != max(idx))]
    res1 <- do.call(cbind, res1)
  }

  vec <- as.vector(res1[9:12,])
  vec <- vec[!is.na(vec)]
  plot(sort(vec), seq(0, 1, length.out = length(vec)), asp = T)
  lines(c(0,1), c(0,1), col = "red")

  power <- length(which(vec <= 0.05/2))/length(vec)
  std <- .bootstrap_power(vec)

  c(power, std)
}


