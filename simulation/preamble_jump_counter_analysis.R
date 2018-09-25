rm(list=ls())
load("../simulation/preamble_jump_counter.RData")
true_jumps <- c(100,140)

bs_trials <- sapply(bs_res, function(x){
  vec <- sapply(1:ncol(x), function(i){
    sum(sapply(1:2, function(k){
      length(which(abs(x[,i] - true_jumps[k]) <= 2))
    }))
  })
  vec <- cumsum(vec)
  which.min(abs(vec - 500))
})
bs_trials # 3800, 2200, 1000, 400, 300, 200

fl_trials <- sapply(fl_res, function(x){
  vec <- sapply(1:ncol(x), function(i){
    sum(sapply(1:2, function(k){
      length(which(abs(x[,i] - true_jumps[k]) <= 2))
    }))
  })
  vec <- cumsum(vec)
  which.min(abs(vec - 500))
})
fl_trials  # 2300, 1300, 600, 400, 300, 200
