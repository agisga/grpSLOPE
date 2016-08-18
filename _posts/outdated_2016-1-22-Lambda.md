---
layout: post
title: Regularizing sequences for Group SLOPE
published: false
---

Currently the package `grpSLOPE` offers the function `lambdaGroupSLOPE` to generate the regularizing sequence $\lambda\subscript{1}, \ldots, \lambda\subscript{J}$ according to one of the following methods (specified via the function argument `method`):

* "BH" &mdash; method of Theorem 1.1 in Bogdan et. al. (2015) *SLOPE &mdash; Adaptive variable selection via convex optimization*.
* "gaussian" &mdash; method of Section 3.2.2 in Bogdan et. al. (2015) *SLOPE &mdash; Adaptive variable selection via convex optimization*.
* "gaussianMC" &mdash; method introduced in [Gossmann et. al. (2015) *Identification of Significant Genetic Variants via SLOPE, and Its Extension to Group SLOPE*](http://dx.doi.org/10.1145/2808719.2808743).
* "chiOrthoMax" &mdash; method of Theorem 2.5 in [Brzyski et. al. (2015) *Group SLOPE — adaptive selection of groups of predictors*](http://arxiv.org/abs/1511.09078).
* "chiOrthoMean" &mdash; $\lambda$s of equation (2.14) in [Brzyski et. al. (2015) *Group SLOPE — adaptive selection of groups of predictors*](http://arxiv.org/abs/1511.09078).
* "chiEqual" &mdash; Procedure 1 in [Brzyski et. al. (2015) *Group SLOPE — adaptive selection of groups of predictors*](http://arxiv.org/abs/1511.09078).
* "chiMean" &mdash; Procedure 2 in [Brzyski et. al. (2015) *Group SLOPE — adaptive selection of groups of predictors*](http://arxiv.org/abs/1511.09078).
*  "chiMC" &mdash; an experimental Monte Carlo based method which is inspired by equation (2.25) [Brzyski et. al. (2015) *Group SLOPE — adaptive selection of groups of predictors*](http://arxiv.org/abs/1511.09078).

For illustration purposes we compute these regularizing sequences for simulated data consisting of 2000 observations with 90 groups of predictors (30 groups of size 5, 30 groups of size 10, and 30 groups of size 20).

The methods "gaussianMC" and "chiMC" require the model matrix $A$. For "gaussianMC" we simulated $A$ from the multivariate normal distribution, such that the mean pairwise correlation of predictors from the same group is 0.9497 (variance: 4.9137e-06), and the average pairwise correlation of variables from different groups is 0.0555 (variance: 0.0067). For "chiMC" we generate a 2000 &times; 1050 random matrix with entries sampled from the normal distribution with mean 0 and variance 1/2000. The columns of this matrix are standardized to have Euclidean norms equal to one before it is used in `lambdaGroupSLOPE`.
Since "chiMC" also requires the response vector $y$, we select 10 significant groups in the coefficient vector $b$ at random and generate $y = Ab + \varepsilon$, where $\varepsilon$ are i.i.d. standard normal variables.

The following code demonstrates how the regularizing sequence can be computed with each of the above methods (except for "chiEqual" because it requires equal group sizes). The first ten coefficients of each of the resulting sequences are graphed right below.

```R
fdr     <- 0.1
n.obs   <- 2000
n.group <- 90
group   <- vector()
for (i in 1:30) {
  tmp <- rep((i-1)*3+c(1,2,3), c(5,10,20))
  group <- c(group, tmp)
}
wt <- rep(c(sqrt(5), sqrt(10), sqrt(20)), 30)
names(wt) <- names(getGroupID(group))

lambda.BH <- lambdaGroupSLOPE(fdr=fdr, n.group=n.group, method="BH")
lambda.G <- lambdaGroupSLOPE(fdr=fdr, n.group=n.group, n.obs=n.obs,
                             method="gaussian")
lambda.MC <- lambdaGroupSLOPE(fdr=fdr, group=group, A=A, n.MC=30, 
                              MC.reps=5000, method="gaussianMC")
lambda.max <- lambdaGroupSLOPE(fdr=fdr, group=group, wt=wt, 
                               method="chiOrthoMax") 
lambda.mean <- lambdaGroupSLOPE(fdr=fdr, group=group, wt=wt, 
                                method="chiOrthoMean") 
lambda.chi <- lambdaGroupSLOPE(fdr=fdr, n.obs=n.obs, group=group, 
                               wt=wt, method="chiMean")
lambda.chiMC <- lambdaGroupSLOPE(fdr=fdr, group=group, A=B, y=y, wt=wt,
                                 n.MC=10, MC.reps=5000, method="chiMC")
```

![Graph of lambda sequences](/grpSLOPE/img/lambda_seq.png)
