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

###################################

rm(list=ls())
load("../experiments/experiment_jump_counter_2.RData")
true_jumps <- c(100,140)

# first a plot of detection
.grab_lines <- function(bs_res, fl_res, paramMat, k = 2){
  idx_vec <- which(paramMat[,"k"] == k)

  bs_line <- sapply(bs_res[idx_vec], function(x){
    sum(apply(x, 2, function(y){
      sum(sapply(true_jumps, function(true_jump){
        length(which(abs(y - true_jump) <= 2))
      }))
    }))/prod(dim(x))
  })

  fl_line <- sapply(fl_res[idx_vec], function(x){
    sum(apply(x, 2, function(y){
      sum(sapply(true_jumps, function(true_jump){
        length(which(abs(y - true_jump) <= 2))
      }))
    }))/prod(dim(x))
  })

  list(bs_line = bs_line, fl_line = fl_line)
}

par(mfrow = c(1,2))
tmp <- .grab_lines(bs_res, fl_res, paramMat, k = 2)
plot(NA, xaxt = "n", xlab='delta', ylab = "cond. powers", xlim = c(1,6),
     ylim = c(0,1), main = "Detection, 2-step algorithms")
axis(1, at=1:6, labels=paramMat[1:6,1])
lines(tmp$bs_line, col = 1, lwd = 2)
lines(tmp$fl_line, col = 3, lwd = 2)

tmp <- .grab_lines(bs_res, fl_res, paramMat, k = 3)
plot(NA, xaxt = "n", xlab='delta', ylab = "cond. powers", xlim = c(1,6),
     ylim = c(0,1), main = "Detection, 3-step algorithms")
axis(1, at=1:6, labels=paramMat[1:6,1])
lines(tmp$bs_line, col = 1, lwd = 2)
lines(tmp$fl_line, col = 3, lwd = 2)

# next a detection numbers
.grab_mat <- function(bs_res, fl_res, paramMat, k = 2, SnR = 2){
  idx <- intersect(which(paramMat[,"k"] == k), which(paramMat[,"SnR"] == SnR))

  bs_vec <- table(apply(bs_res[[idx]], 2, function(x){
    sum(sapply(true_jumps, function(true_jump){
      length(which(abs(x - true_jump) <= 2))
    }))
  }))
  bs_vec <- bs_vec/sum(bs_vec)

  fl_vec <- table(apply(fl_res[[idx]], 2, function(x){
    sum(sapply(true_jumps, function(true_jump){
      length(which(abs(x - true_jump) <= 2))
    }))
  }))
  fl_vec <- fl_vec/sum(fl_vec)

  list(len.fl = fl_vec, len.bsfs = bs_vec)
}

.grab_mat(bs_res, fl_res, paramMat, k = 2, SnR = 2)
.grab_mat(bs_res, fl_res, paramMat, k = 2, SnR = 4)
######
.grab_mat(bs_res, fl_res, paramMat, k = 3, SnR = 2)
.grab_mat(bs_res, fl_res, paramMat, k = 3, SnR = 4)
######
.grab_mat(bs_res, fl_res, paramMat, k = 4, SnR = 2)
.grab_mat(bs_res, fl_res, paramMat, k = 4, SnR = 4)

###############
# last the same example


