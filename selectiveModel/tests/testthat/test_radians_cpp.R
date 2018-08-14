context("Test radians C++")

.intersect_intervals <- function(lis){
  Reduce(.intersect_two_intervals, lis)
}

.intersect_two_intervals <- function(mat1, mat2){
  vec <- sort(unique(c(as.numeric(mat1), as.numeric(mat2))))
  vec <- sort(c(vec, .interpolate(vec)))

  bool_vec <- sapply(vec, function(x){
    all(.theta_in_interval(x, mat1), .theta_in_interval(x, mat2))
  })

  idx_mat <- .consecutive_true(bool_vec)

  mat <- matrix(vec[idx_mat], ncol = 2)
  mat <- mat[order(mat[,1]),,drop = F]

  stopifnot(all(as.numeric(t(mat)) == sort(as.numeric(t(mat)))))
  mat
}

.interpolate <- function(vec){
  n <- length(vec)
  sapply(2:n, function(x){
    mean(vec[(x-1):x])
  })
}

.theta_in_interval <- function(theta, mat){
  stopifnot(abs(theta) <= pi/2)

  vec <- any(apply(mat, 1, function(x){
    x[1] <= theta & theta <= x[2]
  }))
}

.consecutive_true <- function(vec){
  idx <- which(vec)
  if(length(idx) == 0) stop("No intersection")
  breakpoint <- which(sapply(2:length(idx), function(x){idx[x]-idx[x-1] != 1}))

  breakpoint <- c(0, breakpoint, length(idx))
  mat <- t(sapply(2:length(breakpoint), function(x){
    c(idx[breakpoint[x-1]+1], idx[breakpoint[x]])
  }))

  #remove singletons
  idx <- which(mat[,1] != mat[,2])
  if(length(idx) == 0) stop("No intersection")
  mat[idx,]
}

.generate_random_list <- function(i){
  set.seed(10*i)
  num_elements <- ceiling(runif(1, 1, 200))
  common_interval <- sort(runif(2, -pi/2, pi/2))
  mid_point <- mean(common_interval)
  lapply(1:num_elements, function(x){
    .interval(range(c(runif(2, -pi/2 + 0.1, pi/2), common_interval)),
              mid_point)
  })
}

test_that("intersect_intervals_simult is correct", {
  trials <- 100
  bool <- sapply(1:trials, function(x){
    lis <- .generate_random_list(x)

    res1 <- .intersect_intervals(lis)
    res2 <- intersect_intervals_simult(lis)

    all.equal(res1, res2)
  })

  expect_true(all(bool))
})
