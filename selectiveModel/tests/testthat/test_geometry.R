context("Test geometry")

## .quadratic is correct

test_that(".quadratic works", {
  res <- .quadratic(-1, 2, 5)

  expect_true(length(res) == 2)
  expect_true(res[1] < res[2])
  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(all(!is.na(res)))
})

test_that(".quadratic can return only one number", {
  res <- .quadratic(1, 2, 1)

  expect_true(length(res) == 1)
  expect_true(!is.na(res))
})

test_that(".quadratic can return only NA", {
  res <- .quadratic(1, 2, 5)

  expect_true(length(res) == 1)
  expect_true(is.na(res))
})

#################################

## .point_on_plane is correct

test_that(".point_on_plane works", {
  obj <- .plane(c(1,2,3,4), 4)
  res <- .point_on_plane(obj)

  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == 4)
})

test_that(".point_on_plane returns a point on the plane", {
  trials <- 100
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    obj <- .plane(rnorm(10), rnorm(1))
    res <- .point_on_plane(obj)

    ifelse(abs(res%*%obj$a - obj$b) < 1e-6, TRUE, FALSE)
  })

  expect_true(all(bool_vec))
})

test_that(".point_on_plane handles cases parallel to axes", {
  trials <- 100
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    vec <- rep(0, 10)
    idx <- sample(1:10, 5)
    vec[idx] <- rnorm(5)
    obj <- .plane(vec, rnorm(1))
    res <- .point_on_plane(obj)

    ifelse(abs(res%*%obj$a - obj$b) < 1e-6, TRUE, FALSE)
  })

  expect_true(all(bool_vec))
})

test_that(".point_on_plane handles cases where parallel to only but one axes", {
  obj <- .plane(c(0,0,1), 0)
  res <- .point_on_plane(obj)

  expect_true(abs(res%*%obj$a - obj$b) < 1e-6)
})
