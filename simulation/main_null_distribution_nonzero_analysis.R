rm(list=ls())
load("../simulation/main_null_distribution_nonzero.RData")

.plot_qq <- function(res){
  pvalue_vec <- as.numeric(res[[1]][9:12,])
  pvalue_vec <- pvalue_vec[!is.na(pvalue_vec)]
  print(length(pvalue_vec))
  plot(sort(pvalue_vec), seq(0, 1, length.out = length(pvalue_vec)), asp = T)
  lines(c(0,1), c(0,1), col = "red", lwd = 2)

  invisible()
}

.plot_qq(bs_res)
.plot_qq(fl_res)


