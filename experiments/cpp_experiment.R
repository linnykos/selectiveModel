lis = list(matrix(1:10, ncol = 2), matrix(c(1,5,2,10), ncol = 2), matrix(c(1,4,7,10), ncol = 2))
#unlist_native(lis)
theta_in_all_matrix(1,lis)
theta_in_all_matrix(3,lis)

theta_in_matrix(1, lis[[2]])
theta_in_matrix(3, lis[[2]])

####

mat1 <- .partition_interval(c(-7*pi/6, -pi/3))
mat2 <- .partition_interval(c(-pi/2, 0))
mat3 <- .partition_interval(c(-3*pi/4, -pi/4))

res <- .intersect_two_intervals(mat1, mat2)
res <- .intersect_two_intervals(res, mat3)

lis <- list(mat1, mat2, mat3)
# zz = intersect_intervals(lis)
zz = func2()

func2()
