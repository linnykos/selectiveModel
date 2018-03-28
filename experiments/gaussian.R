rm(list=ls())
set.seed(10)
n <- 6
samples <- 1000
lapse <- 50

mean_vec <- c(0,0,0,1,1,1)
cov_mat <- diag(n)
cov_mat[1:3,1:3] <- 1
cov_mat[4:6,4:6] <- 1
diag(cov_mat) <- 2

mat1 <- MASS::mvrnorm(n = samples, mu = mean_vec, Sigma = cov_mat)
mat3 <- MASS::mvrnorm(n = samples, mu = mean_vec, Sigma = cov_mat)


segments <- .segments(6, 1, 1)
gaussian <- .gaussian(mean_vec, cov_mat)
mat2 <- matrix(0, nrow = samples*lapse, ncol = n)
mat2[1,] <- mat1[1,]

for(i in 2:nrow(mat2)){
  #v <- as.numeric(.sample_matrix_space(segments, 1, null = T))
  v <- rnorm(n); v <- v/.l2norm(v)

  rotation <- .rotation_matrix(v, c(1, rep(0, n-1)))
  gaussian2 <- .transform_gaussian(gaussian, mat2[i-1,], rotation)
  univariate <- .conditional_gaussian(gaussian2, rep(0, n-1))
  alpha <- rnorm(1, mean = univariate$mean, sd = sqrt(univariate$covariance))

  mat2[i,] <- mat2[i-1,] + alpha*v
}
mat2 <- mat2[seq(1, nrow(mat2), by = lapse),]

par(mfrow = c(2,3))
for(i in 1:6){
  plot(sort(mat1[,i]), sort(mat2[,i]), asp = T)
  lines(c(-1000,1000), c(-1000,1000), col = "red")
}


par(mfrow = c(2,3))
for(i in 1:6){
  plot(sort(mat1[,i]), sort(mat3[,i]), asp = T)
  lines(c(-1000,1000), c(-1000,1000), col = "red")
}
