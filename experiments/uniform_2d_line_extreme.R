#create a polyhedra that basically includes all points
mat <- rbind(diag(4), -diag(4))
u <- c(rep(-10,4), rep(-10,4))
polyhedra <- structure(list(gamma = mat, u = u), class =)

set.seed(10)
n <- 4
y <- rnorm(n)
y_org <- y
fit <- binSegInf::binSeg_fixedSteps(y, 1)

segments <- .segments(n, binSegInf::jumps(fit))
plane <- .plane(segments, .segment_means(y, segments))
gaussian <- .gaussian(y, diag(n))

tmp <- .sample_matrix_space(segments, null = T)
v <- tmp[,1]; w <- tmp[,2]
