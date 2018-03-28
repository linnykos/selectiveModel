rm(list=ls())
set.seed(10)
n <- 6
trials <- 5000

# determine polyhedron
y <- rnorm(n)
obj <- binSegInf::binSeg_fixedSteps(y, 1)
poly <- binSegInf::polyhedra(obj)

###############

segments <- .segments(length(y), jump_vec = binSegInf::jumps(obj), ignore_jump = 1)
segments_full <- rbind(c(1, rep(0, n-1)), c(0, 1, rep(0, n-2)),
                       c(0,0,1,0,0,0), c(0,0,0,1,0,0), c(0,0,0,0,1,0),
                       segments)
mean_val <- as.numeric(segments%*%y)
cov_mat <- diag(n)

cov_mat_full <- segments_full %*% cov_mat %*% t(segments_full)
mean_cond <- cov_mat_full[1:5,6]%*%solve(cov_mat_full[6,6])%*%mean_val
cov_cond <- cov_mat_full[1:5,1:5] - cov_mat_full[1:5,6]%*%solve(cov_mat_full[6,6])%*%t(cov_mat_full[1:5,6])


.sampler_guassian <- function(samples){
  mat <- matrix(0, ncol = n, nrow = samples)
  i <- 1

  while(i <= samples){
    vec <- MASS::mvrnorm(n = 1, mu = mean_cond, Sigma = cov_cond)
    vec <- solve(segments_full)%*%c(vec, mean_val)

    if(all(poly$gamma%*%vec >= poly$u)){
      mat[i,] <- vec
      i <- i+1
    }
  }

  mat
}

set.seed(10)
samples <- t(.sampler_guassian(trials))
samples3 <- t(.sampler_guassian(trials))

###############
gaussian <- .gaussian(rep(0, n), covariance = diag(n))
set.seed(10)
samples2 <- .sampler_hit_run_line(y, gaussian, segments, poly,
                                 num_samp = trials, burn_in = 100, lapse = 10)


#############

#plot a QQ-plot along each dimension
par(mfrow = c(3,2))
for(i in 1:6){
  plot(sort(samples[i,]), sort(samples2[i,]), asp = T)
  lines(c(-1000,1000), c(-1000,1000), col = "red")
}

par(mfrow = c(3,2))
for(i in 1:6){
  plot(sort(samples3[i,]), sort(samples2[i,]), asp = T)
  lines(c(-1000,1000), c(-1000,1000), col = "red")
}

################

#plot distribution of the test statistic
segments_include <- .segments(length(y), jump_vec = binSegInf::jumps(obj))
dist1 <- apply(segments_include %*% samples, 2, function(x){diff(range(x))})
dist2 <- apply(segments_include %*% samples2, 2, function(x){diff(range(x))})
dist3 <- apply(segments_include %*% samples3, 2, function(x){diff(range(x))})

hist(dist1, breaks = 50, col = "gray")
hist(dist2, breaks = 50, col = rgb(1,0,0,0.5), add = T)
hist(dist3, breaks = 50, col = rgb(0,1,0,0.5), add = T)

par(mfrow = c(1,2))
plot(sort(dist1), sort(dist2), asp = T, pch = 16, col = rgb(0,0,0,0.1))
lines(c(-1000,1000), c(-1000,1000), col = "red")
plot(sort(dist1), sort(dist3), asp = T, pch = 16, col = rgb(0,0,0,0.1))
lines(c(-1000,1000), c(-1000,1000), col = "red")
