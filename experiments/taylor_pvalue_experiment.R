rm(list=ls())
set.seed(10)
source("../experiments/taylor_pvalue.R")

trials <- 2000

pvalue_calculation <- function(seed){
  if(seed %% floor(trials/10) == 0) print('*')
  set.seed(seed)
  n <- 6
  y <- rnorm(n)
  obj <- binSegInf::binSeg_fixedSteps(y, 1)
  poly <- binSegInf::polyhedra(obj)
  segments <- .segments(length(y), jump_vec = binSegInf::jumps(obj), ignore_jump = 1)
  gaussian <- .gaussian(rep(0, n), diag(n))

  samp <- taylor_pvalue_hitrun_sampler(y, gaussian, poly, segments)

  #check
  bool_vec <- apply(samp, 2, function(x){
    all(poly$gamma %*% x >= poly$u)
  })
  stopifnot(all(bool_vec))

  mean_dist <- apply(samp, 2, function(x){
    segments %*% x
  })
  stopifnot(diff(range(mean_dist)) < 1e-6)

  segment_real <- .segments(length(y), jump_vec = binSegInf::jumps(obj))
  segment_dist <- apply(samp, 2, function(x){
    diff(segment_real %*% x)
  })

  length(which(abs(segment_dist) > as.numeric(abs(diff(segment_real %*% y)))))/length(segment_dist)
}

pval_vec <- sapply(1:trials, pvalue_calculation)

plot(sort(pval_vec), seq(0, 1, length.out = trials), asp = T)
lines(c(0,1), c(0,1), col = "red")
