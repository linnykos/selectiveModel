context("Test radian functions")

## .radius is correct

test_that(".radius works", {
  set.seed(10)
  y <- rnorm(10)
  v <- rnorm(10); w <- rnorm(10)
  v <- v/.l2norm(v)
  w <- .projection(w, v); w <- w/.l2norm(w)

  res <- .radius(0, y, v, w)

  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == 10)
})

test_that(".radius is cyclical with period pi, as negatives", {
  trials <- 100
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    y <- rnorm(10)
    v <- rnorm(10); w <- rnorm(10)
    v <- v/.l2norm(v)
    w <- .projection(w, v); w <- w/.l2norm(w)
    theta <- runif(1, min = -pi/2, max = pi/2)

    res1 <- .radius(theta, y, v, w)
    res2 <- .radius(theta+pi, y, v, w)

    ifelse(sum(abs(res1 + res2)) < 1e-6, TRUE, FALSE)
  })

  expect_true(all(bool_vec))
})

################################
