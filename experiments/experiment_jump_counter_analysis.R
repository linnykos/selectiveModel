rm(list=ls())
load("../experiments/experiment_jump_counter.RData")
true_jumps <- c(100,140)

bs_count <- sapply(bs_res, function(x){
  x <- apply(x, 2, sort, decreasing = F)
  sum(sapply(1:2, function(i){length(which(abs(x[i,] - true_jumps[i]) <= 2))}))
})

fl_count <- sapply(fl_res, function(x){
  x <- apply(x, 2, sort, decreasing = F)
  sum(sapply(1:2, function(i){length(which(abs(x[i,] - true_jumps[i]) <= 2))}))
})
