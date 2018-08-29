rm(list=ls())
set.seed(30)
y <- rnorm(10)
obj <- binseginf::binSeg_fixedSteps(y, 2)
poly <- binseginf::polyhedra(obj)
segments <- .segments(length(y), jump_vec = binseginf::jumps(obj))
norm1 <- .l2norm(y)

polyhedra = poly
num_samp = 25
burn_in = 10
lapse = 2
verbose = F

num_col <- burn_in + num_samp*lapse
y_mat <- matrix(NA, nrow = length(y), ncol = num_samp)
seq_idx <- burn_in + (1:num_samp)*lapse
polyhedra <- .remove_rows(polyhedra, .l2norm(y))

prev_y <- y

for(i in 1:35){
  print(i)
  next_y <- .hit_run_next_point_radial(prev_y, segments, polyhedra)
  if(i %in% seq_idx){
    y_mat[,which(seq_idx == i)] <- next_y
  }

  prev_y <- next_y
}

i = 36
y = prev_y
tmp <- .sample_matrix_space(segments, 2, null = T)
v <- tmp[,1]; w <- tmp[,2]

### BRANCH INTO R-CODE
res_list1 <- sapply(1:length(poly$u), function(x){
  plane <- .plane(poly$gamma[x,], poly$u[x])
  plane <- .intersect_plane_basis(plane, y, v, w)
  if(any(is.na(plane))) {
    print(x)
    return(matrix(c(-pi/2, pi/2), ncol = 2))
  }
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

### BRANCH INTO C-CODE
interval_list <- lapply(1:nrow(polyhedra$gamma), function(x){
  print(x)
  c_form_interval(polyhedra$gamma[x,], polyhedra$u[x], y, v, w)
})
