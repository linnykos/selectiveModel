set.seed(6)
n <- 3
y <- rnorm(n)
fit_method <- function(y){binSegInf::binSeg_fixedSteps(y, 1)}

set.seed(1)
res <- selected_model_inference(y, fit_method,
                                test_func = selectiveModel::segment_difference,
                                verbose = F, cores = NA,
                                param = list(burn_in = 1, lapse = 1, time_limit = 600),
                                num_samp = 3000, sigma = 1, ignore_jump = 1)

#plot along 2-dimensional plane
fit <- binSegInf::binSeg_fixedSteps(y, 1)
segments <- .segments(n, binSegInf::jumps(fit), ignore_jump = 1)
v <- .sample_matrix_space(segments, null = T)
v <- cbind(v, as.numeric(segments))

zz <- apply(res$samples, 2, function(x){
  solve(v, x)
})

col_vec <- sapply(1:ncol(res$samples), function(x){
  z <- ncol(res$samples)
  rgb(x/z, 0, (z-x)/z)
})

plot(zz[1,], zz[2,], col = col_vec, pch = 16, asp = T)

xx <- solve(v, y); points(xx[1], xx[2], col = "black", pch = 16, cex = 4)
