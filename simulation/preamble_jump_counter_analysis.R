rm(list=ls())
load("../simulation/preamble_jump_counter_onejump.RData")
true_jumps <- c(160)
target_trials <- 1000

bs_trials <- sapply(bs_res, function(x){
  if(!is.matrix(x)) x <- matrix(x, nrow = 1)

  vec <- sapply(1:ncol(x), function(i){
    sum(sapply(1:2, function(k){
      length(which(abs(x[,i] - true_jumps[k]) <= 2))
    }))
  })
  vec <- cumsum(vec)
  which.min(abs(vec - target_trials))
})
names(bs_trials) <- names(bs_res)
bs_trials # 7600, 4400, 1900, 800, 500, 400


fl_trials <- sapply(fl_res, function(x){
  if(!is.matrix(x)) x <- matrix(x, nrow = 1)

  vec <- sapply(1:ncol(x), function(i){
    sum(sapply(1:2, function(k){
      length(which(abs(x[,i] - true_jumps[k]) <= 0))
    }))
  })
  vec <- cumsum(vec)
  which.min(abs(vec - target_trials))
})
names(fl_trials) <- names(fl_res)
fl_trials  # 4600, 2500, 1250, 650, 450, 350
