set.seed(10)
n <- 4
y <- rnorm(n)
y_org <- y
fit <- binSegInf::binSeg_fixedSteps(y, 1)
polyhedra <- binSegInf::polyhedra(fit)

segments <- .segments(n, binSegInf::jumps(fit))
plane <- .plane(segments, .segment_means(y, segments))

###########

trials <- 1000
theta_vec <- rep(NA, trials)
tmp <- .sample_matrix_space(segments, 2, null = T)
v <- tmp[,1]; w <- tmp[,2]
y_mat <- matrix(NA, ncol = trials, nrow = length(y))
y_mat[,1] <- y

for(trial in 1:(trials-1)){
  interval <- .range_theta_polyhedra(y_mat[,trial], v, w, polyhedra)

  if(nrow(interval) == 1) {
    theta <- stats::runif(1, interval[1], interval[2])
  } else {
    len <- interval[,2] - interval[,1]
    row_idx <- sample(1:nrow(interval), 1, prob = len)
    theta <- stats::runif(1, interval[row_idx,1], interval[row_idx,2])
  }

  theta_vec[trial] <- theta

  y_mat[,trial+1] <- .radians_to_data(theta, y_mat[,trial], v, w)
}

#x_vec <- sin(theta_vec); y_vec <- cos(theta_vec)
#plot(x_vec, y_vec, asp = T)

apply(y_mat, 2, .l2norm)
apply(y_mat, 2, function(x){.segment_means(x, segments)})

y_mat <- y_mat[,-1]
coef_mat <- apply(y_mat, 2, function(x){
  vec <- x - y_org
  mat <- cbind(v, w)
  a <- mat[c(1,3), c(1:2)]
  b <- vec[c(1,3)]
  coef <- solve(a, b)

  #theta <- atan(b/a)
  #if(a < 0) theta <- theta + pi
})

plot(coef_mat[1,], coef_mat[2,], asp = T)
##WARNING: the theta's keep changing so... we need to actually compute the theta wrt the
## original y
