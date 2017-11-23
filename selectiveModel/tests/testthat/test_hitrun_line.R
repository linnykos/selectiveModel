context("Test hit and run, Line")

## .hit_run_next_point_line is correct

test_that(".hit_run_next_point_line works", {
  set.seed(15)
  y <- rnorm(10)
  obj <- binSegInf::binSeg_fixedSteps(y, 2)
  poly <- binSegInf::polyhedra(obj)
  segments <- .segments(length(y), jump_vec = binSegInf::jumps(obj))
  gaussian <- .gaussian(rep(0,10), diag(10))

  res <- .hit_run_next_point_line(y, segments, poly, gaussian)

  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == length(y))
})

test_that(".hit_run_next_point_line preserves the mean", {
  set.seed(5)
  y <- rnorm(10)
  obj <- binSegInf::binSeg_fixedSteps(y, 2)
  poly <- binSegInf::polyhedra(obj)
  segments <- .segments(length(y), jump_vec = binSegInf::jumps(obj))
  gaussian <- .gaussian(rep(0,10), diag(10))

  trials <- 100
  bool_vec <- sapply(1:trials, function(x){
    res <- .hit_run_next_point_line(y, segments, poly, gaussian)

    sum(abs(.segment_means(y, segments) - .segment_means(res, segments))) < 1e-6
  })

  expect_true(all(bool_vec))
})

test_that(".hit_run_next_point_line returns points inside the polyhedra", {
  set.seed(5)
  y <- rnorm(10)
  obj <- binSegInf::binSeg_fixedSteps(y, 2)
  poly <- binSegInf::polyhedra(obj)
  segments <- .segments(length(y), jump_vec = binSegInf::jumps(obj))
  gaussian <- .gaussian(rep(0,10), diag(10))

  trials <- 100
  bool_vec <- sapply(1:trials, function(x){
    res <- .hit_run_next_point_line(y, segments, poly, gaussian)

    all(poly$gamma %*% res >= poly$u)
  })

  expect_true(all(bool_vec))
})

