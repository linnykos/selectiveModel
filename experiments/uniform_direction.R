set.seed(10)

trials <- 2000
v_mat2 <- matrix(NA, ncol = trials, nrow = 3)
poly_gamma <- matrix(c(-1, 0, 0,
                       1, 0, 0,
                       0, -1, 0,
                       0, 1, 0), ncol = 3, nrow = 4, byrow = T)
poly_u <- c(-1, .5, -3, -3)

.l2norm <- function(x){sqrt(sum(x^2))}

basis1 <- c(1,1,1); basis1 <- basis1/.l2norm(basis1)
basis2 <- c(0,-1,1); basis2 <- basis2/.l2norm(basis2)
basis3 <- c(-2,1,1); basis3 <- basis3/.l2norm(basis3)
basis <- rbind(basis1, basis2, basis3)

.sample_from_constraint <- function(basis){
  cov_mat2 <- basis%*%t(basis)
  z <- MASS::mvrnorm(1, mu = rep(0, 2), Sigma = cov_mat2[2:3,2:3]-cov_mat2[2:3,1,drop = F]%*%cov_mat2[1,2:3]/cov_mat2[1,1])
  x <- solve(basis)%*%c(0, z)
}

i <- 1
while(TRUE){
  #generate two samples, both having mean 0
  x1 <- .sample_from_constraint(basis)
  x2 <- .sample_from_constraint(basis)

  #check to see if both samples lie to polyhedron
  if(any(poly_gamma %*% x1 < poly_u) | any(poly_gamma %*% x2 < poly_u)) next()

  #compute the difference, record the direction
  dif <- x1 - x2
  v_mat2[,i] <- dif/.l2norm(dif)

  i <- i + 1
  if(i > trials) break()
}

v_mat2 <- apply(v_mat2, 2, function(x){
  solve(t(basis), x)[2:3]
})

#plot the difference with respect to the null space
pdf("../experiments/uniform_direction.pdf", height = 4, width = 4)
plot(NA, xlim = c(-1, 1), ylim = c(-1, 1), asp = T,
     xlab = "Basis 2", ylab = "Basis 3", main = "Using rejection sampling")
points(v_mat2[1,], v_mat2[2,], pch = 16, col = rgb(0,0,0, 0.075))
graphics.off()
