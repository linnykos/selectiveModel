set.seed(10)

trials <- 2000
v_mat <- matrix(NA, ncol = trials, nrow = 3)
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

start_y <- c(0.75,0,0)
gaussian <- .gaussian(rep(0, 3), covariance = diag(3))
segments <- .segments(3, 1, 1)
polyhedra <- structure(list(gamma = poly_gamma, u = poly_u), class = "polyhedra")

samples <- .sampler_hit_run_line(start_y, gaussian, segments, polyhedra,
                                 num_samp = 2*trials, burn_in = 1000, lapse = 10)

v_mat <- sapply(1:trials, function(x){
  dif <- samples[,x]-samples[,x+trials]
  dif/.l2norm(dif)
})


v_mat <- apply(v_mat, 2, function(x){
  solve(t(basis), x)[2:3]
})

#plot the difference with respect to the null space
pdf("../experiments/uniform_direction2.pdf", height = 4, width = 4)
plot(NA, xlim = c(-1, 1), ylim = c(-1, 1), asp = T,
     xlab = "Basis 2", ylab = "Basis 3", main = "Using Hit-and-run")
points(v_mat[1,], v_mat[2,], pch = 16, col = rgb(0,0,0, 0.075))
graphics.off()

##########

angle_vec <- apply(v_mat, 2, function(x){
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
angle2_vec <- apply(v_mat2, 2, function(x){
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

pdf("../experiments/uniform_direction_qq.pdf", height = 4, width = 4)
plot(sort(angle_vec), sort(angle2_vec), asp = T, pch = 16, col = rgb(0, 0, 0, 0.075),
     xlab = "CDF of angles using rejection sampling", ylab = "CDF of angles using hit-and-run")
lines(c(-100,100),c(-100,100), lwd = 2, col = "red")
graphics.off()
