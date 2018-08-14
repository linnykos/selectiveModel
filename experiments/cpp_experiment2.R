set.seed(10)
n <- 10
y <- c(rep(0, n/2), rep(1, n/2)) + rnorm(n)
fit <- binSegInf::binSeg_fixedSteps(y, 1)
polyhedra <- binSegInf::polyhedra(fit)
gaussian <- .gaussian(rep(0, n), diag(n))

segments <- .segments(n, binSegInf::jumps(fit), 1)

start_y = y;
nullspace_mat <- .sample_matrix_space(segments)
mean_val <- as.numeric(segments%*%start_y)
segments_full <- rbind(t(nullspace_mat), segments)

setting_1 <- .remove_nullspace(gaussian, polyhedra, segments_full, mean_val)
setting_2 <- .whiten(setting_1$gaussian, setting_1$polyhedra)
new_polyhedra <- setting_2$polyhedra

start_z <- setting_2$forward_translation(setting_1$forward_translation(start_y))
n <- length(start_z)
directions <- .generate_directions(n)
alphas <- -new_polyhedra$gamma %*% t(directions)

#flip the slack to accommodate our setting
slack <- -new_polyhedra$gamma %*% start_z + new_polyhedra$u

###################

sample_truncnorm_white(start_z, slack, t(directions), alphas, 10, 10)

state = start_z
direction = t(directions)[,1]
U = slack
alpha = alphas[,1]

state
U

gibbs_step(state, direction, U, alpha)
state
U
