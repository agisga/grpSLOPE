library(grpSLOPE)

context("lambdaGroupSLOPE()")

n.obs <- 50
n.group <- 10
group <- c(10, 4, 7,10, 6, 9, 9, 2, 5, 1, 6, 2, 6, 1, 2, 1,
           2, 1, 5, 9, 5, 1, 3,10, 4, 5, 3, 7, 6, 9)
A <- as.matrix(read.table("./test_data/gaussianMC_test_mat.txt"))
fdr <- 0.1

test_that("method 'BH' computes the sequence correctly", {
  # sol is obtained via the function create_lambda in package SLOPE
  sol <- c(2.575829,2.326348,2.170090,2.053749,1.959964,
           1.880794,1.811911,1.750686,1.695398,1.6448540)
  lambda.BH <- lambdaGroupSLOPE(fdr=fdr, n.group=n.group, method="BH")
  expect_equal(lambda.BH, sol, tolerance=1e-6)
})

test_that("method 'gaussian' computes the sequence correctly", {
  # sol is obtained via the function create_lambda in package SLOPE
  sol <- c(2.575829,2.481928,2.447715,2.437303,2.437303,
           2.437303,2.437303,2.437303,2.437303,2.437303)
  lambda.G <- lambdaGroupSLOPE(fdr=fdr, n.group=n.group, n.obs=n.obs, method="gaussian")
  expect_equal(lambda.G, sol, tolerance=1e-6)
})

test_that("method 'gaussianMC' produces a non-increasing sequence", {
  lambda.MC <- lambdaGroupSLOPE(fdr=fdr, group=group, A=A, n.obs=n.obs,
                                method="gaussianMC", MC.reps=50)
  # check for NA and NaN
  expect_true(!anyNA(lambda.MC))
  # check non-increasing
  expect_true(all(diff(lambda.MC) <= 0))
})
