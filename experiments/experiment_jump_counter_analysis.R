rm(list=ls())
load("../experiments/experiment_jump_counter.RData")
true_jumps <- c(100,140)

bs_count <- sapply(bs_res, function(x){
  x <- apply(x, 2, sort, decreasing = F)
  sum(sapply(1:2, function(i){length(which(abs(x[i,] - true_jumps[i]) <= 2))}))
})
bs_count

fl_count <- sapply(fl_res, function(x){
  x <- apply(x, 2, sort, decreasing = F)
  sum(sapply(1:2, function(i){length(which(abs(x[i,] - true_jumps[i]) <= 2))}))
})
fl_count

##############

bs_trials <- sapply(bs_res, function(x){
  x <- apply(x, 2, sort, decreasing = F)
  seq_vec <- c(1:120)*100
  count <- sapply(seq_vec, function(k){
    sum(sapply(1:2, function(i){length(which(abs(x[i,1:k] - true_jumps[i]) <= 2))}))
  })
  seq_vec[which.min(abs(count - 2000))]
})
bs_trials # 24000, 6600, 2200, 1300, 1100

fl_trials <- sapply(fl_res, function(x){
  x <- apply(x, 2, sort, decreasing = F)
  seq_vec <- c(1:120)*100
  count <- sapply(seq_vec, function(k){
    sum(sapply(1:2, function(i){length(which(abs(x[i,1:k] - true_jumps[i]) <= 2))}))
  })
  seq_vec[which.min(abs(count - 2000))]
})
fl_trials # 13000, 6100, 2900, 1800, 1300

