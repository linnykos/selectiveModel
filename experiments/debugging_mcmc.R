rm(list=ls())
set.seed(44)
n <- 6
trials <- 1000
paramMat <- matrix(c(0,1), nrow = 1)
fit_method <- function(x){binSegInf::binSeg_fixedSteps(x, numSteps = 1)}
vec <- paramMat[1,]
lev <- vec[1]
mean_vec <- c(rep(0, floor(n/2)), rep(lev, floor(n/2)))
dat  <- mean_vec + stats::rnorm(length(mean_vec))
fit <- binSegInf::binSeg_fixedSteps(dat, 1)
poly <- binSegInf::polyhedra(fit)
contrast <- binSegInf::contrast_vector(fit, vec[2])
test_func <- selectiveModel::segment_difference
num_samp <- 2000
cores <- NA
verbose <- T
# selected <- selected_model_inference(dat, fit_method = fit_method, test_func = test_func,
#                                      num_samp = num_samp, ignore_jump = vec[2], cores = cores, verbose = F, param = list(burn_in = 2000, lapse = 1))

ignore_jump <- vec[2]
param <- list(burn_in = 2000, lapse = 1)
y <- dat
n <- length(y)
fit <- fit_method(y)
polyhedra <- binSegInf::polyhedra(fit)
test_stat <- test_func(y, fit, jump = ignore_jump)
segments <- .segments(n, binSegInf::jumps(fit), ignore_jump = ignore_jump)
param <- .fill_in_arguments(param)
# samples <- .sampler_hit_run_radial(y, segments, polyhedra, num_samp = num_samp,
#                                    burn_in = param$burn_in, lapse = param$lapse,
#                                    verbose = verbose)

burn_in <- param$burn_in
lapse <- param$lapse
num_col <- burn_in + num_samp*lapse
y_mat <- matrix(NA, nrow = length(y), ncol = num_samp)
seq_idx <- burn_in + (1:num_samp)*lapse
prev_y <- y
for(i in 1:num_col){
  if(i == 2169) break()
  next_y <- .hit_run_next_point_radial(prev_y, segments, polyhedra)
  if(i %in% seq_idx){
    y_mat[,which(seq_idx == i)] <- next_y
  }
  prev_y <- next_y
}

y_org <- y
y <- prev_y
tmp <- .sample_matrix_space(segments, 2, null = T)
v <- tmp[,1]; w <- tmp[,2]
# interval <- .range_theta_polyhedra(y, v, w, polyhedra)
# if(nrow(interval) == 1) {
#   theta <- stats::runif(1, interval[1], interval[2])
# } else {
#   len <- interval[,2] - interval[,1]
#   row_idx <- sample(1:nrow(interval), 1, prob = len)
#   theta <- stats::runif(1, interval[row_idx,1], interval[row_idx,2])
# }
# y_new <- .radians_to_data(theta, y, v, w)

interval_list <- lapply(1:nrow(polyhedra$gamma), function(x){
  print(x)
  plane <- .plane(polyhedra$gamma[x,], polyhedra$u[x])
  plane <- .intersect_plane_basis(plane, y, v, w)
  if(any(is.na(plane))) return(matrix(c(-pi/2, pi/2), ncol = 2))
  center <- c(-y%*%v, -y%*%w)
  radius <- sqrt(sum(center^2))
  circle <- .circle(center, radius)
  dis <- .distance_point_to_plane(center, plane)

  if(dis >= radius){
    matrix(c(-pi/2, pi/2), ncol = 2)
  } else {
    mat <- .intersect_circle_line(plane, circle)
    stopifnot(nrow(mat) == 2)
    vec <- apply(mat, 2, .euclidean_to_radian, circle = circle)
    init_theta <- .initial_theta(y, v, w)
    .interval(vec, init_theta)
  }
})

x <- 7
plane <- .plane(polyhedra$gamma[x,], polyhedra$u[x])
plane2 <- .intersect_plane_basis(plane, y, v, w)
# if(any(is.na(plane))) return(matrix(c(-pi/2, pi/2), ncol = 2))
center <- c(-y%*%v, -y%*%w)
radius <- sqrt(sum(center^2))
circle <- .circle(center, radius)
dis <- .distance_point_to_plane(center, plane2)
# if(dis >= radius){
#   matrix(c(-pi/2, pi/2), ncol = 2)
# } else {
  mat <- .intersect_circle_line(plane2, circle)
  stopifnot(nrow(mat) == 2)
  vec <- apply(mat, 2, .euclidean_to_radian, circle = circle)
  init_theta <- .initial_theta(y, v, w)
  .interval(vec, init_theta)
#}

# #####checking .intersect_plane_basis
# zz <- 1000
# bool_vec <- sapply(1:zz, function(x){
#   set.seed(x)
#   test_z <- 5*rnorm(2)
#   bool1 <- all(plane2$a %*% test_z >= plane2$b)
#
#   test_y <- y + test_z[1]*v + test_z[2]*w
#   bool2 <- all(plane$a%*%test_y >= plane$b)
#
#   bool1 == bool2
# })
# all(bool_vec)
#
# ######checking .distance_point_to_plane
# zz <- 1000
# bool_vec <- sapply(1:zz, function(x){
#   set.seed(x) #fails for 732
#   test_z <- c(5*rnorm(1), 0)
#   test_z[2] <- (plane2$b - test_z[1]*plane2$a[1])/plane2$a[2]
#   dis2 <- .l2norm(test_z - center)
#
#   dis2 >= dis
# })
# all(bool_vec)

######.intersect_circle_line
bool_vec <- apply(mat, 2, function(x){ #check to be on plane. WARNING: this should be bound by columns, not rows...
  abs(plane2$a %*%x - plane2$b) < 1e-16
})
all(bool_vec)

bool_vec <- apply(mat, 2, function(x){ #check to be on plane
  abs(.l2norm(x - circle$center) - circle$radius) < 1e-16
})
all(bool_vec)

##### bug seems to be in .intersect_circle_line
plane <- plane2
tol <- 1e-6
tol2 <- 1e-9
dis <- .distance_point_to_plane(circle$center, plane)
#if(dis > circle$radius + tol) return(NA)
#if(abs(plane$a[1]) > tol){
#  a1 <- plane$a[1]; a2 <- plane$a[2]
#} else {
  a1 <- plane$a[2]; a2 <- plane$a[1]
#}

a <- 1 + (a2/a1)^2
b <- -2*(a2/a1)*(plane$b/a1 - circle$center[1]) -
  2*circle$center[2]
c <- -circle$radius^2 + (plane$b/a1 - circle$center[1])^2 +
  circle$center[2]^2

y <- .quadratic(a, b, c)
stopifnot(all(!is.na(y)))
x <- (plane$b -a2*y)/a1
