library(Matrix)

n.obs <- 700
n.significant.blocks <- 10 #number of significant blocks
signal <- 1

# (1) Generate the n.obs x 1050 design matrix X from the multivariate normal distribution
# with mean 0 and 1050x1050 covariance matrix Sigma,
# where Sigma is block diagonal with blocks Gamma of sizes
# 5x5, 10x10, and 20x20 (30 blocks of each size),
# such that Gamma[i,i] = 1 if i=j and =0.9 otherwise.
blockcovmat5 <- diag(1,5,5)
blockcovmat5[upper.tri(blockcovmat5,diag=F)] <- blockcovmat5[lower.tri(blockcovmat5, diag=F)] <- 0.95
blockcovmat10 <- diag(1,10,10)
blockcovmat10[upper.tri(blockcovmat10,diag=F)] <- blockcovmat10[lower.tri(blockcovmat10, diag=F)] <- 0.95
blockcovmat20 <- diag(1,20,20)
blockcovmat20[upper.tri(blockcovmat20,diag=F)] <- blockcovmat20[lower.tri(blockcovmat20, diag=F)] <- 0.95
Sigma <- bdiag(rep(list(blockcovmat5, blockcovmat10, blockcovmat20), 30))
withinblock.ind <- as.matrix(Sigma) # for (2) below
withinblock.ind[upper.tri(withinblock.ind,diag=T)] <- 0
#increase between group cor
Sigma <- as.matrix(Sigma)
Sigma[which(Sigma==0)] <- 0.05
#---
Sigma.chol <- as.matrix(chol(Sigma))
A <- matrix(rnorm(mean=0, sd=1/1050, 1050*n.obs),n.obs, 1050) %*% Sigma.chol
# A needs colnorms 1
A <- apply(A, 2, function(x) x/sqrt(sum(x^2)) )

# (2) Compute the average within and between blocks correlations
corA <- abs(cor(A))
withinblock.cor <- list("mean"=mean(c(corA[which(withinblock.ind!=0)])),
                        "var"=var(c(corA[which(withinblock.ind!=0)])))
betweenblock.cor <- list("mean"=mean(c(corA[which(withinblock.ind==0)])),
                         "var"=var(c(corA[which(withinblock.ind==0)])))

# (3) Create the vector b consisting of 90 blocks of sizes 5, 10, and 20
b <- rep(0, 1050)
# indices of the first entries of each block
block.start <- rep(NA, 90)
block.count <- 1
for (i in 1:30){
  block.start[3*i-2] <- block.count
  block.start[3*i-1] <- block.count + 5
  block.start[3*i] <- block.count + 15
  block.count <- block.count + 35
}
# choose significant blocks
nsignificant <- n.significant.blocks
block.significant <- sort(sample(90,nsignificant))
for(i in block.significant){
  if (i == 90){
    b[block.start[90]:1050] <- signal
  }
  else{
    b[block.start[i]:(block.start[i+1]-1)] <- (-1)^i * signal
  }
}

# (4) Create the response vector y
errorvector <- rnorm(n.obs,0,1)
y <- A %*% b + errorvector

# (5) Create B and z ~ B %*% b
B <- matrix(rnorm(700*1050, sd=1/700), 700, 1050)
B.colnorms <- sqrt(apply(B^2, 2, sum))
B <- B %*% diag(1/B.colnorms)
z <- B %*% b + errorvector


# ----------------- Compute lambda sequences --------------------

fdr <- 0.1
n.group <- 90
group <- vector()
for (i in 1:30) {
  tmp <- rep((i-1)*3+c(1,2,3), c(5,10,20))
  group <- c(group, tmp)
}
wt <- rep(c(sqrt(5), sqrt(10), sqrt(20)), 30)
names(wt) <- names(getGroupID(group))

lambda.BH <- lambdaGroupSLOPE(fdr=fdr, n.group=n.group, method="BH")
lambda.G <- lambdaGroupSLOPE(fdr=fdr, n.group=n.group, n.obs=n.obs, method="gaussian")
lambda.MC <- lambdaGroupSLOPE(fdr=fdr, group=group, A=A, n.MC=30, MC.reps=5000, method="gaussianMC")
lambda.max <- lambdaGroupSLOPE(fdr=fdr, group=group, wt=wt, method="chiOrthoMax")
lambda.mean <- lambdaGroupSLOPE(fdr=fdr, group=group, wt=wt, method="chiOrthoMean")
lambda.chi <- lambdaGroupSLOPE(fdr=fdr, n.obs=n.obs, group=group, wt=wt, method="chiMean")
lambda.chiMC <- lambdaGroupSLOPE(fdr=fdr, group=group, A=B, y=z, wt=wt, n.MC=10, MC.reps=5000, method="chiMC")

plot(lambda.BH[1:10], type="l", col=1, ylab="lambda", xlab="index",
     main="Lambda sequence for Group SLOPE",
     ylim=range(c(lambda.BH, lambda.G, lambda.MC, lambda.max, lambda.mean, lambda.chi)))
lines(lambda.G[1:10], col=2)
lines(lambda.MC[1:10], col=3)
lines(lambda.max[1:10], col=4)
lines(lambda.mean[1:10], col=5)
lines(lambda.chi[1:10], col=6)
lines(lambda.chiMC[1:10], col=8)
legend(c("BH", "gaussian", "gaussianMC", "chiOrthoMax", "chiOrthoMean", "chiMean", "chiMC"),
       x=7, y=3.3, col=c(1:6,8), lty=rep(1,7))
