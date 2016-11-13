# grpSLOPE

[![Build Status](https://travis-ci.org/agisga/grpSLOPE.svg?branch=master)](https://travis-ci.org/agisga/grpSLOPE)
[![CRAN downloads this month](http://cranlogs.r-pkg.org/badges/grpSLOPE)](https://CRAN.R-project.org/package=grpSLOPE)

Group SLOPE is a penalized linear regression method that is used for adaptive selection of groups of significant predictors in a high-dimensional linear model. A unique feature of the Group SLOPE method is that it offers (group) false discovery rate control (i.e., control of the expected proportion of irrelevant groups among the total number of groups of predictors selected by Group SLOPE).
A detailed description of the method can be found in [D. Brzyski, A. Gossmann, W. Su, and M. Bogdan (2016) *Group SLOPE &mdash; adaptive selection of groups of predictors*](https://arxiv.org/abs/1610.04960).

## Usage

The basic model fitting procedure for a Group SLOPE model is:

```R
library(grpSLOPE)

result <- grpSLOPE(X=X, y=y, group=group, fdr=0.1)
```

where `X` is the model matrix, `y` the response vector, `group` a vector specifying group memberships of the predictor variables, and `fdr` the target (group) false discovery rate for the variable selection procedure.

One can access various estimated parameters of the fitted model, such as the selected groups, the estimated noise level, or the estimated vector of regression coefficients:

```R
# groups that have a non-zero effect
result$selected
# noise level
sigma(result)
# coefficient vector
coef(result)
```

[A detailed basic usage example can be found here](http://www.alexejgossmann.com/grpSLOPE/basic-usage/).

More complicated (and possibly less helpful) example codes are available in the repository [grpSLOPE_examples](https://github.com/agisga/grpSLOPE_examples).

## Installation

The latest stable version of `grpSLOPE` can be installed from CRAN (The Comprehensive R Archive Network). Just open an R session and do:

```R
install.packages("grpSLOPE")
```

## Installation of the development version

The easiest way to install the latest development version of `grpSLOPE` is by using the R package `devtools`. Just open up an R session and run:

```R
# Install devtools, if you haven't already.
install.packages("devtools")

library(devtools)
install_github("agisga/grpSLOPE")
```

## Contributing to the R package

### Code style

Variable names are all lower case with words separated by dots.
Function names begin with a lower case letter and are written in camel case.
Constants names are all caps.
Otherwise, the package code follows [Google's R style guide](https://google.github.io/styleguide/Rguide.xml).

### Workflow

0. Modify the code.
1. Open [`grpSLOPE.Rproj`](https://github.com/agisga/grpSLOPE/blob/master/grpSLOPE.Rproj) with RStudio.
2. Run `devtools::document()`.
3. Do "Build and Reload" from the menu (or CTRL-Shift-B).
4. Do `devtools::test()` to run the unit tests.
5. Install with `devtools::install()`
6. Run checks with `devtools::check()`
