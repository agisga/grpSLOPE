library(grpSLOPE)

context("proxGroupSortedL1()")

test_that("the prox is evaluated correctly if the groups are consequtive blocks", {
  x   <- proxGroupSortedL1(y = 1:10, group = c(1,1,1,2,2,2,3,3,3,4), lambda = 10:1)
  sol <- c(0, 0, 0, 0.353261553, 0.441576942, 0.529892330, 
           1.974292890, 2.256334731, 2.538376572, 1.0000000001)
  expect_equal(x, sol)
})

test_that("the prox is evaluated correctly if the groups are not consequtive blocks", {
  y   <- 1:10
  grp <- c(1,1,1,2,2,2,3,3,3,4)
  sol <- c(0, 0, 0, 0.353261553, 0.441576942, 0.529892330, 
           1.974292890, 2.256334731, 2.538376572, 1.0000000001)
  ord <- sample(y, 10)
  x <- proxGroupSortedL1(y = y[ord], group = grp[ord], lambda = 10:1)
  expect_equal(x, sol[ord], tolerance=1e-6)
})



context("grpSLOPE()")

test_that("the Group SLOPE solution is computed correctly", {
  set.seed(1)
  A    <- matrix(runif(100, 0, 1), 10, 10)
  grp  <- c(0, 0, 1, 1, 2, 2, 2, 2, 2, 3)
  Dinv <- diag(c(0.5, 0.5, 0.5, 0.5, 0.2, 0.2, 0.2, 0.2, 0.2, 1))
  x    <- c(0, 0, 5, 1, 0, 0, 0, 1, 0, 3)
  y    <- A %*% x
  result <- grpSLOPE(X=A, Dinv=Dinv, b=y, group=grp, lambda=0.1*(10:1),
                   list("tolerance" = 1e-12, "verbosity" = 0))
  solution  <- matrix(c(0,0,3.856002988,2.080742942,0,0,0,0,0,3.512828045), c(10,1))
  expect_equal(result$x, solution, tolerance=1e-6)
})
