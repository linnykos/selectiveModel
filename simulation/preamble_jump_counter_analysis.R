rm(list=ls())
load("../simulation/preamble_jump_counter.RData")
true_jumps <- c(100,140)
target_trials <- 1000

bs_trials <- sapply(bs_res, function(x){
  vec <- sapply(1:ncol(x), function(i){
    sum(sapply(1:2, function(k){
      length(which(abs(x[,i] - true_jumps[k]) <= 2))
    }))
  })
  vec <- cumsum(vec)
  which.min(abs(vec - target_trials))
})
bs_trials # 7600, 4400, 1900, 800, 500, 400

fl_trials <- sapply(fl_res, function(x){
  vec <- sapply(1:ncol(x), function(i){
    sum(sapply(1:2, function(k){
      length(which(abs(x[,i] - true_jumps[k]) <= 0))
    }))
  })
  vec <- cumsum(vec)
  which.min(abs(vec - target_trials))
})
fl_trials  # 4600, 2500, 1250, 650, 450, 350

############################

rm(list=ls())
load("../simulation/preamble_jump_counter.RData")
true_jumps <- c(100,140)
target_trials <- 500

bs_trials <- sapply(bs_res, function(x){
  vec <- sapply(1:ncol(x), function(i){
    if(all(c(100, 140) %in% x[,i])) 1 else 0
  })
  vec <- cumsum(vec)
  which.min(abs(vec - target_trials))
})
bs_trials # 600

fl_trials <- sapply(fl_res, function(x){
  vec <- sapply(1:ncol(x), function(i){
    if(all(c(100, 140) %in% x[,i])) 1 else 0
  })
  vec <- cumsum(vec)
  which.min(abs(vec - target_trials))
})
fl_trials  # 550


