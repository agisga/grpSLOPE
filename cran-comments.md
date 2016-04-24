## Test environments
* local ubuntu 15.10, R 3.2.5
* local BunsenLabs GNU/Linux 8.4, R devel (2016-04-22)
* local OS X 10.9.5, R 3.1.1
* ubuntu 12.04 (on travis-ci), R 3.2.4
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs. 

There were 2 NOTEs:

* This NOTE appears only with linux systems (i.e., it does not appear with OS X or win-builder):
  
  checking installed package size ... NOTE
  installed size is  6.9Mb
  sub-directories of 1Mb or more:
    libs   6.7Mb

  The size is due to a .so file resulting from the usage of RcppEigen.

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Alexej Gossmann <alexej.go@googlemail.com>'
  New submission

  It's my first CRAN submission.

## Resubmission

This is a resubmission. In this version I have changed the title in DESCRIPTION as was suggested by Uwe Ligges.
