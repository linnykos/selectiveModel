context("Test rejection sampling")

## .sampler_rejection is correct

test_that(".sampler_rejection works", {
  set.seed(100)
  y <- rnorm(5)
  obj <- binSegInf::binSeg_fixedSteps(y, 1)
  poly <- binSegInf::polyhedra(obj)
  segments <- .segments(length(y), jump_vec = binSegInf::jumps(obj))

  res <- .sampler_rejection(y, segments, poly, num_samp = 10)

  expect_true(is.numeric(res))
  expect_true(is.matrix(res))
  expect_true(all(dim(res) == c(5,10)))
})

test_that(".sampler_rejection preserves segment means", {
  set.seed(50)
  y <- rnorm(5)
  obj <- binSegInf::binSeg_fixedSteps(y, 1)
  poly <- binSegInf::polyhedra(obj)
  segments <- .segments(length(y), jump_vec = binSegInf::jumps(obj))
  seg_mean <- .segment_means(y, segments)

  res <- .sampler_rejection(y, segments, poly, num_samp = 25)

  bool_vec <- sapply(1:ncol(res), function(x){
    seg_mean2 <- .segment_means(res[,x], segments)
    sum(abs(seg_mean - seg_mean)) < 1e-6
  })

  expect_true(all(bool_vec))
})

test_that(".sampler_rejection preserves l2norm", {
  set.seed(100)
  y <- rnorm(5)
  obj <- binSegInf::binSeg_fixedSteps(y, 1)
  poly <- binSegInf::polyhedra(obj)
  segments <- .segments(length(y), jump_vec = binSegInf::jumps(obj))

  res <- .sampler_rejection(y, segments, poly, num_samp = 25)

  bool_vec <- sapply(1:ncol(res), function(x){
    abs(.l2norm(res[,x]) - .l2norm(y)) < 1e-6
  })

  expect_true(all(bool_vec))
})

test_that(".sampler_rejection are all in polyhedra", {
  set.seed(120)
  y <- rnorm(5)
  obj <- binSegInf::binSeg_fixedSteps(y, 1)
  poly <- binSegInf::polyhedra(obj)
  segments <- .segments(length(y), jump_vec = binSegInf::jumps(obj))

  res <- .sampler_rejection(y, segments, poly, num_samp = 25)

  bool_vec <- sapply(1:ncol(res), function(x){
    all(poly$gamma %*% res[,x] >= poly$u)
  })

  expect_true(all(bool_vec))
})
