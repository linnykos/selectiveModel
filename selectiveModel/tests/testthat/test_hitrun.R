context("Test hit and run")

## .projection is correct

test_that(".projection works", {
  set.seed(10)
  vec1 <- rnorm(10); vec2 <- rnorm(10)
  res <- .projection(vec1, vec2)

  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == 10)
})

test_that(".projection preserves orthogonal vectors", {
  set.seed(1)
  vec1 <- rnorm(10); vec2 <- rnorm(10)
  res <- .projection(vec1, vec2)
  res2 <- .projection(res, vec2)

  expect_true(sum(abs(res - res2)) < 1e-6)
})

test_that(".projection removes parallel vectors", {
  set.seed(1)
  vec1 <- rnorm(10); vec2 <- 2*vec1
  res <- .projection(vec1, vec2)

  expect_true(sum(abs(res)) < 1e-6)
})

test_that(".projection reduces the norm", {
  set.seed(1)
  vec1 <- rnorm(10); vec2 <- rnorm(10)
  res <- .projection(vec1, vec2)

  expect_true(.l2norm(vec1) >= .l2norm(res))
})

################################

## .projection_matrix is correct

test_that(".projection_matrix works", {
  set.seed(10)
  mat <- .segments(10, c(3, 7))
  vec <- rnorm(10)
  res <- .projection_matrix(vec, mat)

  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == 10)
})

test_that(".projection_matrix returns a vector that is orthogonal to rows of mat", {
  set.seed(10)
  mat <- .segments(10, c(3, 7))
  vec <- rnorm(10)
  res <- .projection_matrix(vec, mat)

  for(i in 1:nrow(mat)){
    expect_true(abs(res %*% mat[i,]) < 1e-6)
  }
})

test_that(".projection_matrix reduces the norm", {
  set.seed(10)
  mat <- matrix(rnorm(30), ncol = 10, nrow = 3)
  vec <- rnorm(10)
  res <- .projection_matrix(vec, mat)

  expect_true(.l2norm(vec) >= .l2norm(res))
})

test_that(".projection_matrix preserves orthogonal vectors", {
  set.seed(10)
  mat <- matrix(rnorm(30), ncol = 10, nrow = 3)
  vec <- rnorm(10)
  res <- .projection_matrix(vec, mat)
  res2 <- .projection_matrix(res, mat)

  expect_true(sum(abs(res - res2)) < 1e-6)
})

test_that(".projection_matrix removes parallel vectors", {
  set.seed(10)
  mat <- matrix(rnorm(30), ncol = 10, nrow = 3)
  vec <- 2*mat[1,]
  res <- .projection_matrix(vec, mat)

  expect_true(sum(abs(res)) < 1e-6)
})

#########################

## .sample_nullspace is correct

test_that(".sample_nullspace works", {
  set.seed(10)
  mat <- .segments(10, c(3, 7))
  res <- .sample_nullspace(mat, 2)

  expect_true(is.numeric(res))
  expect_true(is.matrix(res))
  expect_true(all(dim(res) == c(10,2)))
})

test_that(".sample_nullspace gives vectors that are orthogonal to mat", {
  set.seed(10)
  mat <- .segments(10, c(3, 7))
  res <- .sample_nullspace(mat, 3)

  for(i in 1:3){
    for(j in 1:3){
      expect_true(abs(res[,i]%*%mat[j,]) < 1e-6)
    }
  }
})

test_that(".sample_nullspace gives vectors that are orthogonal", {
  set.seed(10)
  mat <- .segments(10, c(3, 7))
  vec_mat <- .sample_nullspace(mat, 3)
  res <- t(vec_mat)%*%vec_mat

  expect_true(sum(abs(res - diag(3))) < 1e-6)
})

##############################

## .hit_run_next_point is correct

test_that(".hit_run_next_point works", {
  set.seed(15)
  y <- rnorm(10)
  obj <- binSegInf::binSeg_fixedSteps(y, 2)
  poly <- binSegInf::polyhedra(obj)
  segments <- .segments(length(y), jump_vec = binSegInf::jumps(obj))

  res <- .hit_run_next_point(y, segments, poly)

  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == length(y))
})

test_that(".hit_run_next_point gives a sample with the correct properties", {
  set.seed(10)
  y <- rnorm(10)
  obj <- binSegInf::binSeg_fixedSteps(y, 2)
  poly <- binSegInf::polyhedra(obj)
  segments <- .segments(length(y), jump_vec = binSegInf::jumps(obj))

  res <- .hit_run_next_point(y, segments, poly)

  expect_true(sum(abs(.segment_means(y, segments) -
                        .segment_means(res, segments))) < 1e-6)
  expect_true(abs(.l2norm(y) - .l2norm(res)) < 1e-6)
})

###########################

## .sampler_hit_run is correct

test_that(".sampler_hit_run works", {
  set.seed(50)
  y <- rnorm(10)
  obj <- binSegInf::binSeg_fixedSteps(y, 2)
  poly <- binSegInf::polyhedra(obj)
  segments <- .segments(length(y), jump_vec = binSegInf::jumps(obj))

  res <- .sampler_hit_run(y, segments, poly, num_samp = 25)

  expect_true(is.numeric(res))
  expect_true(is.matrix(res))
  expect_true(all(dim(res) == c(length(y), 25)))
})

test_that(".sampler_hit_run preserves segment means", {
  set.seed(10)
  y <- rnorm(10)
  obj <- binSegInf::binSeg_fixedSteps(y, 2)
  poly <- binSegInf::polyhedra(obj)
  segments <- .segments(length(y), jump_vec = binSegInf::jumps(obj))
  seg_mean <- .segment_means(y, segments)

  res <- .sampler_hit_run(y, segments, poly, num_samp = 25)

  bool_vec <- sapply(1:25, function(x){
    seg_mean2 <- .segment_means(res[,x], segments)
    sum(abs(seg_mean - seg_mean)) < 1e-6
  })

  expect_true(all(bool_vec))
})

test_that(".sampler_hit_run preserves l2norm", {
  set.seed(30)
  y <- rnorm(10)
  obj <- binSegInf::binSeg_fixedSteps(y, 2)
  poly <- binSegInf::polyhedra(obj)
  segments <- .segments(length(y), jump_vec = binSegInf::jumps(obj))
  norm1 <- .l2norm(y)

  res <- .sampler_hit_run(y, segments, poly, num_samp = 25)

  bool_vec <- sapply(1:25, function(x){
    norm2 <- .l2norm(res[,x])
    abs(norm1 - norm2) < 1e-6
  })

  expect_true(all(bool_vec))
})
