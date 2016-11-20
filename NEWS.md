This is a list of API changes.

# grpSLOPE 0.2.0.9000

* Fix for the installation error on r-oldrel-windows-ix86+x86_64 (R v3.2.5) in the CRAN package check results. This error was caused by the generic S3 method `sigma()` not being available from the `stats` package prior to R v3.3.0.

* Checks for missing data in the inputs X, y and group were added in grpSLOPE().

* A check for whether the input matrix X has columns with 0 variance is performed in grpSLOPE() when normalize=TRUE. This should prevent division by 0, when the columns of X are standardized to have norms equal to 1.
