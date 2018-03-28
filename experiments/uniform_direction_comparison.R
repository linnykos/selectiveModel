trials <- 1000

basis1 <- c(1,1,1); basis1 <- basis1/.l2norm(basis1)
basis2 <- c(0,-1,1); basis2 <- basis2/.l2norm(basis2)
basis3 <- c(-2,1,1); basis3 <- basis3/.l2norm(basis3)
basis <- rbind(basis1, basis2, basis3)

segments <- .segments(3, 1, 1)

#method 1
mat1 <- sapply(1:trials, function(x){
  set.seed(x)
  .sample_matrix_space(segments, num_vec = 1)
})

mat1 <- apply(mat1, 2, function(x){
  solve(t(basis), x)[2:3]
})

#method2
mat2 <- sapply(1:trials, function(x){
  set.seed(x)
  vec <- stats::rnorm(3)
  vec <- .projection_matrix(vec, segments)
  vec/.l2norm(vec)
})

mat2 <- apply(mat2, 2, function(x){
  solve(t(basis), x)[2:3]
})

par(mfrow = c(2,1))
plot(mat1[1,], mat1[2,], asp = T, pch = 16, col = rgb(0,0,0,0.05))
plot(mat2[1,], mat2[2,], asp = T, pch = 16, col = rgb(0,0,0,0.05))

###############

angle_vec <- apply(mat1, 2, function(x){
  if(x[1] >= 0 & x[2] >= 0) {
    atan(x[2]/x[1])
  } else if (x[1] <= 0 & x[2] >= 0){
    atan(x[2]/x[1])+pi
  } else if (x[1] <= 0 & x[2] <= 0){
    atan(x[2]/x[1])+pi
  } else {
    atan(x[2]/x[1])+2*pi
  }
})

#now from uniform_direction.R
angle2_vec <- apply(mat2, 2, function(x){
  if(x[1] >= 0 & x[2] >= 0) {
    atan(x[2]/x[1])
  } else if (x[1] <= 0 & x[2] >= 0){
    atan(x[2]/x[1])+pi
  } else if (x[1] <= 0 & x[2] <= 0){
    atan(x[2]/x[1])+pi
  } else {
    atan(x[2]/x[1])+2*pi
  }
})

plot(sort(angle_vec), sort(angle2_vec), asp = T, pch = 16, col = rgb(0, 0, 0, 0.075),
     xlab = "CDF of angles using rejection sampling", ylab = "CDF of angles using hit-and-run")
lines(c(-100,100),c(-100,100), lwd = 2, col = "red")

plot(sort(angle_vec), seq(0,2*pi,length.out = length(angle_vec)),
     asp = T, pch = 16, col = rgb(0, 0, 0, 0.075))
lines(c(-100,100),c(-100,100), lwd = 2, col = "red")
