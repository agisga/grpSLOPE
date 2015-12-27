library(grpSLOPE)

context("Proximal gradient method")

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
  expect_equal(x, sol[ord])
})
