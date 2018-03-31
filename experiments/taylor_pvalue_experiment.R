rm(list=ls())
set.seed(10)
source("../experiments/taylor_pvalue.R")

pvalue_calculation <- function(seed){
  set.seed(seed)
  n <- 6
  y <- rnorm(n)
  obj <- binSegInf::binSeg_fixedSteps(y, 1)
  poly <- binSegInf::polyhedra(obj)
  segments <- .segments(length(y), jump_vec = binSegInf::jumps(obj), ignore_jump = 1)
  gaussian <- .gaussian(rep(0, n), diag(n))

  samp <- taylor_pvalue_hitrun_sampler(y, gaussian, poly, segments)
  segment_real <-  .segments(length(y), jump_vec = binSegInf::jumps(obj))
  segment_dist <- apply(samp, 2, function(x){
    diff(segment_real %*% x)
  })

  length(which(abs(segment_dist) > as.numeric(abs(diff(segment_real %*% y)))))/length(segment_dist)
}

pval_vec <- sapply(1:500, pvalue_calculation)

plot(sort(pval_vec), seq(0, 1, length.out = length(pval_vec)), asp = T)
lines(c(0,1), c(0,1), col = "red")
