library(grpSLOPE)

context("proxGroupSortedL1()")

y   <- 1:10
grp <- c(1,1,1,2,2,2,3,3,3,4)
sol <- c(0, 0, 0, 0.353261553, 0.441576942, 0.529892330, 
         1.974292890, 2.256334731, 2.538376572, 1.0000000001)
lam <- 10:1

test_that("the prox is evaluated correctly when the groups are consequtive blocks", {
  x   <- proxGroupSortedL1(y = y, group = grp, lambda = lam)
  expect_equal(x, sol)
})

test_that("the prox is evaluated correctly when the groups are not consequtive blocks", {
  ord <- sample(y, 10)
  x <- proxGroupSortedL1(y = y[ord], group = grp[ord], lambda = lam)
  expect_equal(x, sol[ord], tolerance=1e-6)
})



context("grpSLOPE()")

set.seed(1)
A    <- matrix(runif(100, 0, 1), 10, 10)
grp  <- c(0, 0, 1, 1, 2, 2, 2, 2, 2, 3)
wt   <- c(0.5, 0.5, 0.5, 0.5, 0.2, 0.2, 0.2, 0.2, 0.2, 1)
x    <- c(0, 0, 5, 1, 0, 0, 0, 1, 0, 3)
y    <- A %*% x
lam  <- 0.1 * (10:1)
sol  <- c(0,0,3.856002988,2.080742942,0,0,0,0,0,3.512828045)

test_that("the Group SLOPE solution is computed correctly when the groups are consequtive blocks", {
  result <- grpSLOPESolver(X=A, Dinv=diag(wt), b=y, group=grp, lambda=lam,
                   list("tolerance" = 1e-12, "verbosity" = 0))
  expect_equal(result$x, as.matrix(sol), tolerance=1e-6)
})

test_that("the Group SLOPE solution is computed correctly when the groups are not consequtive blocks", {
  ord <- sample(1:10, 10)
  result <- grpSLOPESolver(X=A[ , ord], Dinv=diag(wt[ord]), b=y, group=grp[ord], lambda=lam,
                   list("tolerance" = 1e-12, "verbosity" = 0))
  expect_equal(result$x, as.matrix(sol[ord]), tolerance=1e-6)
})
