rm(list=ls())
load("../simulation/main_null_distribution.RData")

.plot_qq <- function(res){
  pvalue_vec <- as.numeric(res[[1]][5:8,])
  pvalue_vec <- pvalue_vec[!is.na(pvalue_vec)]
  plot(seq(0, 1, length.out = length(pvalue_vec)), sort(pvalue_vec), asp = T)
  lines(c(0,1), c(0,1), col = "red", lwd = 2)

  invisible()
}

.plot_qq(bs_res)
.plot_qq(fl_res)


