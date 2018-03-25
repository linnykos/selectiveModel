set.seed(10)

trials <- 1000
v_mat <- matrix(NA, ncol = trials, nrow = 2)
poly_gamma <- matrix(c(1, 0,
                       -1, 0,
                       0, 1,
                       0, -1), ncol = 2, nrow = 4, byrow = T)
poly_u <- c(1, -.5, 3, 3)

.l2norm <- function(x){sqrt(sum(x^2))}

i <- 1
while(TRUE){
  #generate two samples, both having mean 0
  x1 <- rnorm(2)
  x2 <- rnorm(2)

  #check to see if both samples lie to polyhedron
  if(any(poly_gamma %*% x1 > poly_u) | any(poly_gamma %*% x2 > poly_u)) next()

  #compute the difference, record the direction
  dif <- x1 - x2
  v_mat[,i] <- dif/.l2norm(dif)

  i <- i + 1
  if(i > trials) break()
}

#plot the difference with respect to the null space
plot(NA, xlim = c(-1, 1), ylim = c(-1, 1), asp = T)
points(v_mat[1,], v_mat[2,], pch = 16, col = rgb(0.5, 0.5, 0.5, 0.5),
       xlab = "Basis 2", ylab = "Basis 3")
