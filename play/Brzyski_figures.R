# Figure 1 ----------------------

fdr <- 0.1
n.iter <- 300
p <- 5000
X <- diag(rep(1,p))
n.group <- 1000
group <- c(rep(1:200, each=3),
           rep(201:400, each=4),
           rep(401:600, each=5),
           rep(601:800, each=6),
           rep(801:1000, each=7))
group.id <- getGroupID(group)
group.length <- sapply(group.id, FUN=length)
wt <- rep(NA, p)
for (j in 1:n.group) {
  wt[group.id[[j]]] <- sqrt(group.length[j])
}
#wt <- sqrt(group.length)

Bfun <- function(l) {
  sqrt(4*log(n.group) * (1 - n.group^(-2/l)) - l)
}
a <- sum(Bfun(group.length)) / sum(sqrt(group.length))

n.relevant <- floor(seq(1, 250, length=11))

for (k in 1:length(n.relevant)) {
  for (i in 1:n.iter) {
    # generate coeffient vector, pick relevant groups at random
    b <- rep(0, p)
    ind.relevant <- sample(1:n.group, n.relevant[k])
    for (j in ind.relevant) { b[group.id[[j]]] <-a }

    # generate the response vector
    y <- X %*% b + rnorm(p, sd=1)

    # generate lambda
    lambda.BH <- lambdaGroupSLOPE(fdr=fdr, n.group=n.group, method="BH")

    # get Group SLOPE solution
    b.grpSLOPE <- proximalGradientSolverGroupSLOPE(y=y, A=X, group=group, wt=wt, lambda=lambda.BH)

    # FDR and power
    nonzero <- rep(NA, n.group)
    for (j in 1:n.group) { nonzero[j] <- (sum(b.grpSLOPE$x[group.id[[j]]]^2) > 0) }
    truepos <- sum(nonzero[ind.relevant])
    falsepos <- sum(nonzero) - truepos
    FDR <- falsepos / max(1, sum(nonzero))
    pow <- truepos / length(ind.relevant)
    print(paste("FDR:", FDR))
    print(paste("Power:", pow))
  }
}
    
