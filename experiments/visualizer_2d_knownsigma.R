### 3-dimensional example. run binseg, fix the mean, so all the samples lie on a 2d-plane.???

set.seed(6)
n <- 3
y <- rnorm(n)
fit_method <- function(y){binSegInf::binSeg_fixedSteps(y, 1)}

res <- selected_model_inference(y, fit_method,
                                test_func = selectiveModel::segment_difference,
                                verbose = F, cores = NA,
                                num_samp = 250, sigma = 1, ignore_jump = 1,
                                param = list(time_limit = 600))

########################

ymax <- max(c(y, as.numeric(res$samples)))
ymin <- min(c(y, as.numeric(res$samples)))

plot(y, pch = 16, ylim = c(ymin, ymax))
for(i in 1:250){
  points(res$samples[,i], pch = 16, col = rgb(0.5, 0.5, 0.5, 0.1))
  lines(res$samples[,i], pch = 16, col = rgb(0.5, 0.5, 0.5, 0.1))
}
points(y, pch = 16, cex = 2)
points(res$samples[,i], col = "red", pch = 16)


#########################

#plot along 2-dimensional plane
fit <- binSegInf::binSeg_fixedSteps(y, 1)
segments <- .segments(n, binSegInf::jumps(fit), ignore_jump = 1)
v <- .sample_matrix_space(segments, null = T)
v <- cbind(v, as.numeric(segments))

zz <- apply(res$samples, 2, function(x){
  solve(v, x)
})

plot(zz[1,], zz[2,], col = rgb(0.5, 0.5, 0.5, 0.1), pch = 16, asp = T)

xx <- solve(v, y); points(xx[1], xx[2], col = "red", pch = 16)
