rm(list=ls())
#set.seed(100)
vec <- runif(3)
vec <- vec/.l2norm(vec)
A <- matrix(vec, nrow = 1, ncol = 3)

trials <- 1000
vec_mat <- sapply(1:trials, function(x){
  #set.seed(x)
  .sample_matrix_space(A, 1, null = T)
})

vec_mat <- t(vec_mat)
basis1 <- vec_mat[1,]/.l2norm(vec_mat[1,])
basis2 <- .projection(vec_mat[2,], basis1)
basis2 <- basis2/.l2norm(basis2)
basis <- cbind(basis1, basis2)
basis_expanded <- cbind(basis, vec)

stopifnot(abs(basis1 %*% basis2) < 1e-6)

coef <- solve(basis_expanded, t(vec_mat))

plot(coef[1,], coef[2,], asp = T, pch = 16, col = rgb(0.5, 0.5, 0.5, 0.5))
