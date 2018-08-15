context("Test hit and run, Radial")

## .hit_run_next_point_radial is correct

test_that(".hit_run_next_point_radial works", {
  set.seed(15)
  y <- rnorm(10)
  obj <- binSegInf::binSeg_fixedSteps(y, 2)
  poly <- binSegInf::polyhedra(obj)
  segments <- .segments(length(y), jump_vec = binSegInf::jumps(obj))

  res <- .hit_run_next_point_radial(y, segments, poly)

  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == length(y))
})

test_that(".hit_run_next_point_radial gives a sample with the correct properties", {
  set.seed(10)
  y <- rnorm(10)
  obj <- binSegInf::binSeg_fixedSteps(y, 2)
  poly <- binSegInf::polyhedra(obj)
  segments <- .segments(length(y), jump_vec = binSegInf::jumps(obj))

  res <- .hit_run_next_point_radial(y, segments, poly)

  expect_true(sum(abs(.segment_means(y, segments) -
                        .segment_means(res, segments))) < 1e-6)
  expect_true(abs(.l2norm(y) - .l2norm(res)) < 1e-6)
})

###########################

## .sampler_hit_run_radial is correct

test_that(".sampler_hit_run_radial works", {
  set.seed(50)
  y <- rnorm(10)
  obj <- binSegInf::binSeg_fixedSteps(y, 2)
  poly <- binSegInf::polyhedra(obj)
  segments <- .segments(length(y), jump_vec = binSegInf::jumps(obj))

  res <- .sampler_hit_run_radial(y, segments, poly, num_samp = 25, burn_in = 10)

  expect_true(is.numeric(res))
  expect_true(is.matrix(res))
  expect_true(all(dim(res) == c(length(y), 25)))
})

test_that(".sampler_hit_run_radial preserves segment means", {
  set.seed(10)
  y <- rnorm(10)
  obj <- binSegInf::binSeg_fixedSteps(y, 2)
  poly <- binSegInf::polyhedra(obj)
  segments <- .segments(length(y), jump_vec = binSegInf::jumps(obj))
  seg_mean <- .segment_means(y, segments)

  res <- .sampler_hit_run_radial(y, segments, poly, num_samp = 25, burn_in = 10)

  bool_vec <- sapply(1:25, function(x){
    seg_mean2 <- .segment_means(res[,x], segments)
    sum(abs(seg_mean - seg_mean)) < 1e-6
  })

  expect_true(all(bool_vec))
})

test_that(".sampler_hit_run_radial preserves l2norm", {
  set.seed(30)
  y <- rnorm(10)
  obj <- binSegInf::binSeg_fixedSteps(y, 2)
  poly <- binSegInf::polyhedra(obj)
  segments <- .segments(length(y), jump_vec = binSegInf::jumps(obj))
  norm1 <- .l2norm(y)

  res <- .sampler_hit_run_radial(y, segments, poly, num_samp = 25, burn_in = 10)

  bool_vec <- sapply(1:25, function(x){
    norm2 <- .l2norm(res[,x])
    abs(norm1 - norm2) < 1e-6
  })

  expect_true(all(bool_vec))
})

test_that(".sampler_hit_run_radial are all in polyhedra", {
  set.seed(10)
  y <- rnorm(10)
  obj <- binSegInf::binSeg_fixedSteps(y, 2)
  poly <- binSegInf::polyhedra(obj)
  segments <- .segments(length(y), jump_vec = binSegInf::jumps(obj))
  seg_mean <- .segment_means(y, segments)

  res <- .sampler_hit_run_radial(y, segments, poly, num_samp = 25, burn_in = 10)

  bool_vec <- sapply(1:25, function(x){
    .try_polyhedra(res[,x], poly)
  })

  expect_true(all(bool_vec))
})

#####################

## .remove_rows is correct

test_that(".remove_rows works", {
  set.seed(15)
  y <- rnorm(10)
  obj <- binSegInf::binSeg_fixedSteps(y, 2)
  poly <- binSegInf::polyhedra(obj)

  res <- .remove_rows(poly, .l2norm(y))

  expect_true(class(res) == "polyhedra")
  expect_true(ncol(res$gamma) == ncol(poly$gamma))
  expect_true(length(res$u) == nrow(res$gamma))
  expect_true(length(res$u) <= length(poly$gamma))
})

test_that(".remove_rows actually removes rows", {
  mat <- diag(1/1:10)
  vec <- rep(1, 10)
  poly <- binSegInf::polyhedra(mat, vec)

  res <- .remove_rows(poly, 5)

  expect_true(ncol(res$gamma) == ncol(poly$gamma))
  expect_true(length(res$u) < length(poly$gamma))
})

test_that(".remove_rows removes a subset of rows when the radius shrinks", {
  set.seed(15)
  y <- rnorm(10)
  obj <- binSegInf::binSeg_fixedSteps(y, 2)
  poly <- binSegInf::polyhedra(obj)
  poly$u <- rnorm(length(poly$u))

  res1 <- .remove_rows(poly, .l2norm(y))
  res2 <- .remove_rows(poly, .l2norm(y)/10)

  expect_true(nrow(res1$gamma) >= nrow(res2$gamma))
  expect_true(all(res2$u %in% res1$u))
})

test_that(".remove_rows with radius of 0 removes all rows", {
  set.seed(15)
  y <- rnorm(10)
  obj <- binSegInf::binSeg_fixedSteps(y, 2)
  poly <- binSegInf::polyhedra(obj)

  res <- .remove_rows(poly, 0)
  len <- length(which(abs(poly$u) < 1e-6))

  expect_true(nrow(res$gamma) == len)
  expect_true(length(res$u) == len)
})
