library(grpSLOPE)

context("orthogonalizeGroups()")

A <- as.matrix(read.table("./test_data/gaussianMC_test_mat.txt"))

# add an additional group of size 1 to A
A <- cbind(A, 1:nrow(A))

group  <- c(10, 4, 7,10, 6, 9, 9, 2, 5, 1, 6, 2, 6, 1, 2, 1,
            2, 1, 5, 9, 5, 1, 3,10, 4, 5, 3, 7, 6, 9, 8)
grp.id <- getGroupID(group)

test_that("for each group i returns P, Q, R such that A_i[ , P] = Q %*% R", {
  grp.ortho <- orthogonalizeGroups(A, grp.id)
  for (i in 1:length(grp.id)) {
    Ai <- as.matrix(A[ , grp.id[[i]]])
    P  <- grp.ortho[[i]]$P
    Q  <- grp.ortho[[i]]$Q
    R  <- grp.ortho[[i]]$R
    expect_equal(as.matrix(Ai[ , P]), Q %*% R)
  }
})
