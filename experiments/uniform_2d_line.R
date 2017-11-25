set.seed(10)
n <- 4
y <- rnorm(n)
y_org <- y
fit <- binSegInf::binSeg_fixedSteps(y, 1)
polyhedra <- binSegInf::polyhedra(fit)

segments <- .segments(n, binSegInf::jumps(fit))
plane <- .plane(segments, .segment_means(y, segments))
gaussian <- .gaussian(y, diag(n))

tmp <- .sample_matrix_space(segments, null = T)
v <- tmp[,1]; w <- tmp[,2]

#########
# run sampler burn-in for 500 steps first
burn_in <- 500
y2 <- y

for(j in 1:burn_in){
  y2 <- .hit_run_next_point_line(y2, segments, polyhedra, gaussian)
}

#######
# collect samples for real now
trials <- 1000
y_mat <- matrix(NA, ncol = trials, nrow = length(y))
for(j in 1:trials){
  y2 <- .hit_run_next_point_line(y2, segments, polyhedra, gaussian)
  y_mat[,j] <- y2
}

#####
# check
tmp <- apply(y_mat, 2, function(x){.segment_means(x, segments)})
stopifnot(all(apply(tmp, 1, function(x){diff(range(x))}) < 1e-6))

######
# plot

basis <- cbind(v, w)
coef_mat <- apply(y_mat, 2, function(x){
  vec <- x - y
  a <- basis[c(1,3), c(1:2)]
  b <- vec[c(1,3)]
  coef <- solve(a, b)
})

plot(coef_mat[1,], coef_mat[2,], pch = 16, asp = T)

