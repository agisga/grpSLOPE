# Group SLOPE

[![Build Status](https://travis-ci.org/agisga/grpSLOPE.svg?branch=master)](https://travis-ci.org/agisga/grpSLOPE)

Group SLOPE is a penalized linear regression method that is used for adaptive selection of groups of significant predictors in a high-dimensional linear model. It was introduced in [Brzyski et. al. (2015) *Group SLOPE &mdash; adaptive selection of groups of predictors*](http://arxiv.org/abs/1511.09078) and [Gossmann et. al. (2015) *Identification of Significant Genetic Variants via SLOPE, and Its Extension to Group SLOPE*](http://dx.doi.org/10.1145/2808719.2808743).
A unique feature of the Group SLOPE method is that it offers (group) false discovery rate control (i.e., control of the expected proportion of irrelevant groups among the total number of groups of predictors selected by the Group SLOPE method).

## Installation

### Before installation

Your R configuration must allow for a working Rcpp. This is generally not a problem on Unix/Linux, but setting it up on Windows may require some work.

### Installation with devtools (recommended)

The easiest way to install the latest development version of `grpSLOPE` is by using the R package `devtools`. Just open up an R session and run:

```R
# Install devtools, if you haven't already.
install.packages("devtools")

library(devtools)
install_github("agisga/grpSLOPE")
```

### Installation without devtools (not recommended)

If you don't want to use `devtools`, you can install `grpSLOPE` by downloading the source code and then following these steps:

0. Install the R packages `Rcpp` and `RcppEigen` if you don't have them installed already.
1. Go to the directory that contains the `grpSLOPE` directory (which contains the `grpSLOPE` source code).
2. Open an R session and run `Rcpp::compileAttributes("./grpSLOPE")`. Then quit R.
3. Run `R CMD build grpSLOPE`. You should then have a file like `grpSLOPE_0.1.0.tar.gz`.
4. Run `R CMD INSTALL grpSLOPE_0.1.0.tar.gz` to install the package.

## Usage

The basic model fitting procedure is:

```R
library(grpSLOPE)

result <- grpSLOPE(X=X, y=y, group=group, fdr=0.1)
```

where `X` is the model matrix, `y` the response vector, `group` a vector specifying group memberships of the predictor variables, and `fdr` the target (group) false discovery rate for the variable selection procedure.

A detailed basic usage example can be found [here](http://www.alexejgossmann.com/grpSLOPE/basic-usage/).

More complicated (and less helpful) example codes (of varying quality and readability) are available in the repository [grpSLOPE_examples](https://github.com/agisga/grpSLOPE_examples).

## Contributing

### Code style

Variable names are all lower case with words separated by dots.
Function names begin with a lower case letter and are written in camel case.
Constants names are all caps.
Otherwise, I try to follow [Google's R style guide](https://google.github.io/styleguide).

### Workflow

0. Modify the code.
1. Open [`grpSLOPE.Rproj`](https://github.com/agisga/grpSLOPE/blob/master/grpSLOPE.Rproj) with RStudio.
2. Run `devtools::document()`.
3. Do "Build and Reload" from the menu (or CTRL-Shift-B).
4. Do `devtools::test()` to run the unit tests.
5. Install with `devtools::install()`
