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

    ifelse(abs(obj$a%*%res- obj$b) < 1e-6, TRUE, FALSE)
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

    ifelse(abs(obj$a%*%res - obj$b) < 1e-6, TRUE, FALSE)
  })

  expect_true(all(bool_vec))
})

test_that(".point_on_plane handles cases where parallel to only but one axes", {
  obj <- .plane(c(0,0,1), 0)
  res <- .point_on_plane(obj)

  expect_true(abs(obj$a%*%res - obj$b) < 1e-6)
})

test_that(".point_on_plane can handle intersection of planes", {
  mat <- .segments(10, c(3, 7))
  plane <- .plane(mat, 1:3)

  res <- .point_on_plane(plane)

  expect_true(sum(abs(plane$a %*% res - plane$b)) < 1e-6)
})

test_that(".point_on_plane works if point is supposed to be negative",{
  mat <- .segments(5, 3)
  plane <- .plane(mat, c(-0.1498594, 0.5018780))

  res <- .point_on_plane(plane)

  expect_true(sum(abs(plane$a %*% res - plane$b)) < 1e-6)
})

##########################

## .distance_point_to_plane is correct

test_that(".distance_point_to_plane works", {
  plane <- .plane(c(0,0,1), 0)
  point <- c(1,0,1)
  res <- .distance_point_to_plane(point, plane)

  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == 1)
  expect_true(res >= 0)
})

test_that(".distance_point_to_plane is always non-negative", {
  trials <- 100
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    plane <- .plane(rnorm(10), rnorm(1))
    point <- rnorm(10)
    res <- .distance_point_to_plane(point, plane)

    ifelse(res >= 0, TRUE, FALSE)
  })

  expect_true(all(bool_vec))
})

test_that(".distance_point_to_plane satisifies triangle inequality", {
  trials <- 100
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    plane <- .plane(rnorm(10), rnorm(1))
    point1 <- rnorm(10)
    point2 <- rnorm(10)
    res1 <- .distance_point_to_plane(point1, plane)
    res2 <- .distance_point_to_plane(point2, plane)
    dist <- .l2norm(point1 - point2)

    ifelse(abs(res2-res1) <= dist, TRUE, FALSE)
  })

  expect_true(all(bool_vec))
})

####################

## .intersect_circle_line is correct

test_that(".intersect_circle_line works", {
  plane <- .plane(c(0,1), 0)
  circle <- .circle(c(0,0), 1)
  res <- .intersect_circle_line(plane, circle)

  expect_true(is.matrix(res))
  expect_true(is.numeric(res))
  expect_true(all(dim(res) == c(2,2)))
})

test_that(".intersect_circle_line gives a proper point", {
  trials <- 100
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    plane <- .plane(rnorm(2), b = 0)
    circle <- .circle(c(0,0), 1)
    points <- .intersect_circle_line(plane, circle)

    #check to see points are in circle
    res1 <- apply(points, 2, function(x){sum(abs(sqrt(x[1]^2+x[2]^2) - 1)) < 1e-6})

    #check to see points are on plane
    res2 <- apply(points, 2, function(x){abs(plane$a%*%x) < 1e-6})

    all(c(res1, res2))
  })

  expect_true(all(bool_vec))
})

test_that(".intersect_circle_line gives a proper point with b", {
  trials <- 100
  radius <- 100
  b <- 3
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    plane <- .plane(rnorm(2), b = b)
    circle <- .circle(c(0,0), radius)
    points <- .intersect_circle_line(plane, circle)

    #check to see points are in circle
    res1 <- apply(points, 2, function(x){sum(abs(sqrt(x[1]^2+x[2]^2) - radius)) < 1e-6})

    #check to see points are on plane
    res2 <- apply(points, 2, function(x){abs(plane$a%*%x - plane$b) < 1e-6})

    all(c(res1, res2))
  })

  expect_true(all(bool_vec))
})

test_that(".intersect_circle_line gives a proper point with small values of a[1]", {
  trials <- 100
  radius <- 100
  b <- 2
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    a_vec <- c(0, rnorm(1))
    plane <- .plane(a_vec, b = b)
    circle <- .circle(c(0,0), radius*plane$b)
    points <- .intersect_circle_line(plane, circle)

    if(!is.matrix(points)) points <- as.matrix(points, nrow = 1)

    #check to see points are in circle
    res1 <- apply(points, 2, function(x){sum(abs(sqrt(x[1]^2+x[2]^2) - circle$radius)) < 1e-6})

    #check to see points are on plane
    res2 <- apply(points, 2, function(x){abs(plane$a%*%x - plane$b) < 1e-6})

    all(c(res1, res2))
  })

  expect_true(all(bool_vec))
})

test_that(".intersect_circle_line gives a proper point with small values of a[1] always gives 2 points", {
  trials <- 100
  radius <- 100
  b <- 2
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    a_vec <- c(0, rnorm(1))
    plane <- .plane(a_vec, b = b)
    circle <- .circle(c(0,0), radius*plane$b)
    points <- .intersect_circle_line(plane, circle)

    ncol(points) == 2
  })

  expect_true(all(bool_vec))
})

test_that(".intersect_circle_line gives a proper point with small values of a[1] and offcenter circle", {
  trials <- 100
  radius <- 100
  b <- 2
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    a_vec <- c(0, rnorm(1))
    plane <- .plane(a_vec, b = b)
    circle <- .circle(c(1,2), radius*plane$b)
    points <- .intersect_circle_line(plane, circle)

    if(!is.matrix(points)) points <- as.matrix(points, nrow = 1)

    #check to see points are in circle
    res1 <- apply(points, 2, function(x){sum(abs(sqrt((x[1]-circle$center[1])^2+(x[2]-circle$center[2])^2) - circle$radius)) < 1e-6})

    #check to see points are on plane
    res2 <- apply(points, 2, function(x){abs(plane$a%*%x - plane$b) < 1e-6})

    all(c(res1, res2))
  })

  expect_true(all(bool_vec))
})

test_that(".intersect_circle_line gives a proper point with small values of a[2] and offcenter circle", {
  trials <- 100
  radius <- 100
  b <- 2
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    a_vec <- c(rnorm(1), 0)
    plane <- .plane(a_vec, b = b)
    circle <- .circle(c(1,2), radius*plane$b)
    points <- .intersect_circle_line(plane, circle)

    if(!is.matrix(points)) points <- as.matrix(points, nrow = 1)

    #check to see points are in circle
    res1 <- apply(points, 2, function(x){sum(abs(sqrt((x[1]-circle$center[1])^2+(x[2]-circle$center[2])^2) - circle$radius)) < 1e-6})

    #check to see points are on plane
    res2 <- apply(points, 2, function(x){abs(plane$a%*%x - plane$b) < 1e-6})

    all(c(res1, res2))
  })

  expect_true(all(bool_vec))
})

test_that(".intersect_circle_line can give NA", {
  plane <- .plane(c(1,-1), b = 10)
  circle <- .circle(c(0,0), 1)
  res <- .intersect_circle_line(plane, circle)

  expect_true(length(res) == 1)
  expect_true(is.na(res))
})

test_that(".intersect_circle_line can give one point", {
  plane <- .plane(c(1,0), b = 1)
  circle <- .circle(c(0,0), 1)
  res <- .intersect_circle_line(plane, circle)

  expect_true(is.matrix(res))
  expect_true(is.numeric(res))
  expect_true(all(dim(res) == c(2,1)))
})

######################################

## .intersect_plane_basis is correct

test_that(".intersect_plane_basis works", {
  set.seed(10)
  plane <- .plane(rnorm(10), rnorm(1))
  y <- .point_on_plane(plane)
  v <- rnorm(10); w <- rnorm(10)
  v <- v/.l2norm(v)
  w <- .projection(w, v); w <- w/.l2norm(w)

  res <- .intersect_plane_basis(plane, y, v, w)
  expect_true(class(res) == "plane")
  expect_true(length(res$a) == 2)
  expect_true(length(res$b) == 1)
})

test_that(".intersect_plane_basis has the correct properties", {
  set.seed(5)
  plane <- .plane(rnorm(10), -abs(rnorm(1)))
  y <- rep(0, 10)
  v <- rnorm(10); w <- rnorm(10)
  v <- v/.l2norm(v)
  w <- .projection(w, v); w <- w/.l2norm(w)

  res <- .intersect_plane_basis(plane, y, v, w)

  trials <- 100
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    vec <- rnorm(2)
    bool1 <- all(res$a %*% vec >= res$b)

    y_new <- y + vec[1]*v + vec[2]*w
    bool2 <- all(plane$a%*%y_new >= plane$b)

    bool1 == bool2
  })

  expect_true(all(bool_vec))
})

#################################

## .closest_point_to_origin is correct

test_that(".closest_point_to_origin works", {
  mat <- .segments(10, c(3, 7))
  plane <- .plane(mat, 1:3)

  res <- .closest_point_to_origin(plane)

  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == 10)
})

test_that(".closest_point_to_origin returns a point on the plane", {
  mat <- .segments(10, c(3, 7))
  plane <- .plane(mat, 1:3)

  res <- .closest_point_to_origin(plane)

  expect_true(sum(abs(plane$a%*%res - plane$b)) < 1e-6)
})

test_that(".closest_point_to_origin returns the correct point", {
  mat <- .segments(10, c(2,6,9))
  plane <- .plane(mat, 1:4)

  res <- .closest_point_to_origin(plane)
  dis <- .l2norm(res)

  trials <- 100
  nullspace <- .sample_matrix_space(mat)
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    coef <- rnorm(ncol(nullspace))
    vec <- res + nullspace%*%coef

    new_dis <- .l2norm(vec)
    new_dis > dis
  })

  expect_true(all(bool_vec))
})

test_that(".closest_point_to_origin works for this test case", {
  mat <- .segments(5, 3)
  plane <- .plane(mat, c(-0.1498594, 0.5018780))

  res <- .closest_point_to_origin(plane)

  expect_true(sum(abs(plane$a%*%res - plane$b)) < 1e-6)
})

##############################

## .intersect_polyhedron_line is correct

test_that(".intersect_polyhedron_line works", {
  set.seed(15)
  y <- rnorm(10)
  obj <- binSegInf::binSeg_fixedSteps(y, 1)
  poly <- binSegInf::polyhedra(obj)
  line <- .line(y, rep(1, 10))

  res <- .intersect_polyhedron_line(poly, line)

  expect_true(length(res) == 2)
  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
})

test_that(".intersect_polyhedron_line can return two non-Infinite boundaries", {
  set.seed(15)
  y <- c(rep(0,5), rep(5,5)) + rnorm(10)
  obj <- binSegInf::binSeg_fixedSteps(y, 4)
  poly <- binSegInf::polyhedra(obj)
  line <- .line(y, rnorm(10))

  res <- .intersect_polyhedron_line(poly, line)

  expect_true(all(abs(res) <= 1e5))
})

test_that(".intersection_polyhedron_line returns a correct interval", {
  set.seed(10)
  y <- c(rep(0,5), rep(5,5)) + rnorm(10)
  obj <- binSegInf::binSeg_fixedSteps(y, 4)
  poly <- binSegInf::polyhedra(obj)
  line <- .line(y, rnorm(10))

  res <- .intersect_polyhedron_line(poly, line)

  trials <- 100
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    alpha <- runif(1, res[1], res[2])
    vec <- line$point + alpha*line$direction
    all(poly$gamma %*% vec >= poly$u)
  })

  expect_true(all(bool_vec))
})

test_that(".intersection_polyhedron_line gives shorter intervals as more jumps are added", {
  set.seed(10)
  y <- 1:10 + rnorm(10)
  line <- .line(y, rnorm(10))

  obj1 <- binSegInf::binSeg_fixedSteps(y, 1)
  poly1 <- binSegInf::polyhedra(obj1)
  res1 <- .intersect_polyhedron_line(poly1, line)

  obj2 <- binSegInf::binSeg_fixedSteps(y, 5)
  poly2 <- binSegInf::polyhedra(obj2)
  res2 <- .intersect_polyhedron_line(poly2, line)

  expect_true(res1[1] <= res2[1])
  expect_true(res1[2] >= res2[2])
})

test_that(".intersection_polyhedron_line works on harder case", {
  poly <- list(gamma = matrix(c(1,-1, 2,1, -2/3,1, -1/3,-1), nrow = 4, byrow = T),
               u = c(-1, -2, -2, -1))
  y <- c(2, -1/2); v <- c(1,1)
  line <- .line(y, v)
  interval <- .intersect_polyhedron_line(poly, line)

  c <- y + interval[1]*line$direction
  d <- y + interval[2]*line$direction

  expect_true(all(poly$gamma %*% c + 1e-6 >= poly$u))
  expect_true(all(poly$gamma %*% d + 1e-6 >= poly$u))

  c1 <- y + (interval[1]-0.05)*v
  d1 <- y + (interval[2]+0.05)*v

  expect_true(!all(poly$gamma %*% c >= poly$u))
  expect_true(!all(poly$gamma %*% d >= poly$u))
})

test_that(".intersection_polyhedron_line solution remains the same after rotation", {
  set.seed(5)
  y <- rnorm(10)
  obj <- binSegInf::binSeg_fixedSteps(y, 2)
  polyhedra <- binSegInf::polyhedra(obj)
  segments <- .segments(length(y), jump_vec = binSegInf::jumps(obj))
  gaussian <- .gaussian(rep(0,10), diag(10))

  n <- length(y)
  v <- as.numeric(.sample_matrix_space(segments, 1, null = T))
  stopifnot(abs(.l2norm(v)-1) < 1e-6)
  line <- .line(y, v)
  interval <- .intersect_polyhedron_line(polyhedra, line)

  tol <- 1e-6
  rotation <- .rotation_matrix(v, c(1, rep(0, n-1)))
  expect_true(all(polyhedra$gamma %*% (t(rotation)%*%c(1, rep(0, n-1))*(interval[1]+tol) + y) >=
                    polyhedra$u))
  expect_true(all(polyhedra$gamma %*% (t(rotation)%*%c(1, rep(0, n-1))*(interval[2]-tol) + y) >=
                    polyhedra$u))
})


