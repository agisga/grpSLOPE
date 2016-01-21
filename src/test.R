library(Rcpp)
library(RcppEigen)

source("../R/utils.R")

set.seed(1)
grp.id <- getGroupID(sample(1:4,20,repl=T))
X <- matrix(NA, 25, 20)
for (i in 1:4) {
  X[ , grp.id[[i]]] <- as.double(names(grp.id)[i])
}
X <- X + matrix(rnorm(500, sd=1/25), 25, 20)
beta <- as.double(20:1)
y <- X %*% beta + rnorm(25)
lambda <- c(2.2,1.1)

sourceCpp("lambdaMC.cpp")
lambdaChiMCAdjustment(y=y, X=X, group_id=grp.id, lambda=lambda, number_of_drawings=2)
