# Group SLOPE

[![Travis-CI Build Status](https://travis-ci.com/agisga/grpSLOPE.svg?token=MrDRsseBxwE6bVK4mNQS&branch=master)](https://travis-ci.com/agisga/grpSLOPE)

## Installation

The easiest way to install the latest development version of `grpSLOPE` is by using the R package `devtools`. Just open up an R session and run:

```R
# Install devtools, if you haven't already.
install.packages("devtools")

library(devtools)
install_github("agisga/grpSLOPE")
```

If you don't want to use `devtools`, you can install `grpSLOPE` by downloading the source code and then following these steps:

1. Go to the directory that contains the `grpSLOPE` source code directory.
2. Open an R session and run `Rcpp::compileAttributes("./grpSLOPE")`. Then quit R.
3. Run `R CMD build grpSLOPE`. You should then have a file like `grpSLOPE_0.1.0.tar.gz`.
4. Run `R CMD INSTALL grpSLOPE_0.1.0.tar.gz` to install the package.

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
