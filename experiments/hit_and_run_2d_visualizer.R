.plot_poly <- function(poly){
  plot(NA, xlim = c(-1, 3), ylim = c(-2, 1), asp = T)
  points(0, 0, col = "red", cex = 2, pch = 16)
  lines(x = c(-5, 5), y = rep(0, 2), col = "red")
  lines(y = c(-5, 5), x = rep(0, 2), col = "red")
  for(i in 1:nrow(poly$gamma)){
    x_vec <- c(-4, 4)
    y_vec <- (poly$u[i] - poly$gamma[i,1]*x_vec)/poly$gamma[i,2]
    lines(x_vec, y_vec, lwd = 2)
  }
}

create_prob_grid <- function(xlim, ylim, density_func, grid_size = 50){
  xseq <- seq(xlim[1], xlim[2], length.out = grid_size)
  yseq <- seq(ylim[1], ylim[2], length.out = grid_size)

  mat <- sapply(yseq, function(y){ sapply(xseq, function(x){density_func(x,y)})})

  rownames(mat) <- yseq; colnames(mat) <- xseq
  mat/sum(mat)
}


plot_grid <- function(mat, type = c("contour"), alpha = 0.5, nlevels = 10, ...){
  xseq <- as.numeric(colnames(mat)); yseq <- as.numeric(rownames(mat))

  contour(xseq, yseq, mat, nlevels = nlevels, add = T, drawlabels = F, ...)
}

#####################
# design polyhedron
poly <- list(gamma = matrix(c(1,-1, 2,1, -2/3,1, -1/3,-1), nrow = 4, byrow = T),
             u = c(-1, -2, -2, -1))

.plot_poly(poly)

# gaussian contour
gaussian <- .gaussian(mean = c(1, -1), covariance = matrix(c(2,-1,-1,2),2,2))
density_func <- function(x,y){
  mvtnorm::dmvnorm(c(x,y), gaussian$mean, gaussian$covariance)
}
xlim <- c(-1, 3); ylim = c(-2, 1)
mat <- create_prob_grid(xlim, ylim, density_func)
plot_grid(mat, asp = T)

# let's draw points
y <- c(2, -1/2); v <- c(1,1)
line <- .line(y, v)
interval <- .intersect_polyhedron_line(poly, line)
rotation <- .rotation_matrix(v, c(1, rep(0, 1)))
c <- y + interval[1]*line$direction
d <- y + interval[2]*line$direction
points(c[1], c[2], col = "blue", pch = 16)
points(d[1], d[2], col = "green", pch = 16)

gaussian2 <- .transform_gaussian(gaussian, c, rotation)

univariate <- .conditional_gaussian(gaussian2, rep(0, 1))

newsamp <- matrix(NA, nrow = 1000, ncol = 2)
intermediate <- rep(NA, 1000)
for(i in 1:1000){
  intermediate[i] <- .sampler_truncated_gaussian(univariate, 0, interval[2]-interval[1])
  newsamp[i,] <- y + (intermediate[i]+interval[1])*line$direction
}

# plot them
points(newsamp, pch = 16, col = rgb(0.5, 0.5, 0.5, 0.01))
hist(intermediate)
