library(grpSLOPE)

context("proxSortedL1Rcpp()")

y   <- as.double(12:2)
lam <- log(y)
sol <- c(9.515093350,8.602104727,7.697414907,6.802775423,5.920558458, 5.054089851,
         4.208240531,3.390562088,2.613705639,1.901387711,1.306852819)

test_that("the prox is evaluated correctly", {
  x <- proxSortedL1Rcpp(y = y, lambda = lam)
  expect_equal(x, sol)
})
