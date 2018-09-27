rm(list=ls())
load("../simulation/main_null_distribution_nonzero.RData")
load("../simulation/main_null_distribution_nonzero_tmp.RData")

.plot_qq <- function(res){
  pvalue_vec <- as.numeric(res[9:12,])
  pvalue_vec <- pvalue_vec[!is.na(pvalue_vec)]
  print(length(pvalue_vec))
  plot(sort(pvalue_vec), seq(0, 1, length.out = length(pvalue_vec)), asp = T)
  lines(c(0,1), c(0,1), col = "red", lwd = 2)

  invisible()
}

bs_res <- res
bs_res[[1]] <- do.call(cbind, bs_res[[1]][which(sapply(bs_res[[1]], length) > 1)])
.plot_qq(bs_res[[1]])
.plot_qq(bs_res[[1]])


fl_res[[1]] <- do.call(cbind, fl_res[[1]][which(sapply(fl_res[[1]], length) > 1)])
.plot_qq(fl_res[[1]])

idx <- which(fl_res[[2]][5:8,] %in% c(1, -1))
pvalue_vec <- fl_res[[2]][9:12,][idx]
plot(sort(pvalue_vec), seq(0, 1, length.out = length(pvalue_vec)), asp = T)
lines(c(0,1), c(0,1), col = "red", lwd = 2)
.plot_qq(fl_res[[2]])


