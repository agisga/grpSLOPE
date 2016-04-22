## Test environments
* local ubuntu 15.10, R 3.2.5
* ubuntu 12.04 (on travis-ci), R 3.2.4
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs. 

There were 2 NOTEs:

* checking installed package size ... NOTE
  installed size is  7.7Mb
  sub-directories of 1Mb or more:
    libs   7.6Mb 

  The size is due to a .so file resulting from the usage of Rcpp and RcppEigen.

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Alexej Gossmann <alexej.go@googlemail.com>'
  New submission

  It's my first CRAN submission.
