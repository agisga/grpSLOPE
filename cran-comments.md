## Test environments
* local ubuntu 16.04, R 3.3.2 (2016-10-31)
* local BunsenLabs GNU/Linux 8.4, R Under development (unstable) (2016-10-18 r71535)
* x86_64-pc-linux-gnu (64-bit) using R version 3.2.5 (2016-04-14)
* ubuntu 12.04.5 LTS (on travis-ci), R 3.3.1 (2016-06-21)
* win-builder (devel and release)

## R CMD check results
There were no ERRORs, WARNINGs, or NOTEs. 

## CRAN-related changes
Fix for the installation error on r-oldrel-windows-ix86+x86_64 (R v3.2.5) in the CRAN package check results. This error was caused by the generic S3 method `sigma()` not being available from the `stats` package prior to R v3.3.0.
