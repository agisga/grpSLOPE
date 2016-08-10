---
layout: post
title: Basic usage of the R package grpSLOPE
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Basic usage of grpSLOPE}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

## Introduction

Group SLOPE is a penalized linear regression method that is used for adaptive selection of groups of significant predictors in a high-dimensional linear model. It was introduced in [Brzyski et. al. (2015) *Group SLOPE &mdash; adaptive selection of groups of predictors*](http://arxiv.org/abs/1511.09078) and [Gossmann et. al. (2015) *Identification of Significant Genetic Variants via SLOPE, and Its Extension to Group SLOPE*](http://dx.doi.org/10.1145/2808719.2808743).
A unique feature of the Group SLOPE method is that it offers (group) false discovery rate control (i.e., control of the expected proportion of irrelevant groups among the total number of groups of predictors selected by the Group SLOPE method).

## Data generation

We simulate a SNP-data-like model matrix. 

**Note:** In fact, the Group SLOPE method is designed to work with a data matrix, where the columns corresponding to different groups of predictors are nearly uncorrelated. Only for the brevity of exposition we do not check for or enforce low between group correlations in this example.


{% highlight r %}
set.seed(1)

p     <- 500
probs <- runif(p, 0.1, 0.5)
probs <- t(probs) %x% matrix(1,p,2)
X0    <- matrix(rbinom(2*p*p, 1, probs), p, 2*p)
X     <- X0 %*% (diag(p) %x% matrix(1,2,1))
{% endhighlight %}

Upper left $10 \times 10$ corner of $X$:


{% highlight r %}
pander::pandoc.table(X[1:10, 1:10])
{% endhighlight %}


- - - - - - - - - -
0 1 1 1 1 2 2 0 2 1

0 2 0 1 0 2 0 1 1 1

0 1 1 0 0 2 2 1 0 0

1 0 0 1 1 1 1 0 0 0

0 0 2 1 0 1 0 2 1 0

1 0 1 1 0 0 2 1 2 0

0 0 1 1 1 1 1 1 1 2

1 0 0 1 0 0 1 0 2 1

1 0 1 1 0 1 2 0 0 0

0 0 1 1 1 2 1 0 0 0

- - - - - - - - - -

We divide the predictors into 100 groups of sizes ranging from 3 to 7.


{% highlight r %}
group <- c(rep(1:20, each=3),
           rep(21:40, each=4),
           rep(41:60, each=5),
           rep(61:80, each=6),
           rep(81:100, each=7))
group.id <- grpSLOPE::getGroupID(group)
n.group <- length(group.id)
group.length <- sapply(group.id, FUN=length)
{% endhighlight %}

We randomly select 10 groups to be truly significant.


{% highlight r %}
ind.relevant <- sample(1:n.group, 10)
print(sort(ind.relevant))
{% endhighlight %}



{% highlight text %}
##  [1]  4 19 24 26 28 42 50 55 67 74
{% endhighlight %}

Then we generate the vector of predictor coefficients and the response vector according to a linear model.


{% highlight r %}
b <- rep(0, p)
for (j in ind.relevant) {
  b[group.id[[j]]] <- runif(group.length[j])
}

# generate the response vector
y <- X %*% b + rnorm(p)
{% endhighlight %}

## Fitting the Group SLOPE model

We fit the Group SLOPE model to the simulated data. The function argument `fdr` signifies the target group-wise false discovery rate (gFDR) of the variable selection procedure.

{% highlight r %}
library(grpSLOPE)

result <- grpSLOPE(X=X, y=y, group=group, fdr=0.1)
{% endhighlight %}

## Model fit results

We can display which groups were selected as significant by the Group SLOPE method.


{% highlight r %}
result$selected
{% endhighlight %}

{% highlight text %}
##  [1] "4"  "10" "19" "24" "26" "28" "42" "50" "55" "67" "74"
{% endhighlight %}

Similarly, we can look at the estimates of the noise level and the regression coefficients (note that the estimated parameters are returned on the scale of the normalized versions of `X` and `y`).
  

{% highlight r %}
result$sigma
{% endhighlight %}

{% highlight text %}
## [1] 0.9751941
{% endhighlight %}

{% highlight r %}
result$beta[1:14]
{% endhighlight %}

{% highlight text %}
##  [1] 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
##  [8] 0.000000 0.000000 7.374246 2.030661 1.071715 0.000000 0.000000
{% endhighlight %}

It might also be interesting to plot the first few elements of the regularizing sequence $\lambda$ used by the Group SLOPE method for the given inputs.


{% highlight r %}
plot(result$lambda[1:10], xlab = "Index", ylab = "Lambda", type="l")
{% endhighlight %}

![Plot of the lambda sequence](/grpSLOPE/img/2016-2-12-basic-usage/unnamed-chunk-8-1.png)


We check the performance of the method by computing the resulting group false discovery proportion (gFDP) and power.


{% highlight r %}
n.selected    <- length(result$selected)
true.relevant <- names(group.id)[ind.relevant]
truepos       <- intersect(result$selected, true.relevant)

n.truepos  <- length(truepos)
n.falsepos <- n.selected - n.truepos

FDP <- n.falsepos / max(1, n.selected)
pow <- n.truepos / length(true.relevant)

print(paste("FDP =", FDP))
{% endhighlight %}

[1] "FDP = 0.0909090909090909"


{% highlight r %}
print(paste("Power =", pow))
{% endhighlight %}

[1] "Power = 1"

We see that the method indeed did not exceed the target gFDR, while maintaining a high power.

## Lambda sequence

Multiple ways to select the regularizing sequence $\lambda$ are available.

If a group structure with little correlation between groups can be assumed (i.e., groups in the standardized model matrix are nearly orthogonal), then we suggest to use the sequence `corrected`, which is the default.

The sequences `mean` and `max` can be used together with the options `orthogonalize=FALSE` and `normalize=FALSE`, when the columns of the model matrix are exactly orthogonal to each other.

Alternatively, any non-increasing sequence of appropriate length can be used. However, we do not suggest to use any other sequences unless you are an expert on the (Group) SLOPE method.
