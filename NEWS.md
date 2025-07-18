# grpSLOPE (development version)

This is a list of software changes.

# grpSLOPE 0.3.4

* The dependency on the `SLOPE` package in `Suggests` has been removed (see PR #15) due to some breaking changes in `SLOPE`.

# grpSLOPE 0.3.3

* No API changes
* Very minor changes based on feedback from CRAN review.

# grpSLOPE 0.3.2

* No API changes
* Only changes are to the automated software tests to accommodate minor floating point arithmetic differences on M1mac

# grpSLOPE 0.3.1

* Only changes are to the automated software tests.

# grpSLOPE 0.3.0

* Removes the dependency on the R package `SLOPE`, because the only two functions from `SLOPE` that are used in `grpSLOPE` are being deprecated and removed from `SLOPE` package versions newer than 0.1.3.
* Since the dependency on the R package `SLOPE` has been removed (see above), the two functions `SLOPE::SLOPE_solver` and `SLOPE::prox_sorted_L1` (including the underlying C implementation) have been adapted from `SLOPE` version 0.1.3 into this version of the `grpSLOPE` package.
* In addition to the default FISTA solver (function `proximalGradientSolverGroupSLOPE` used by default within function `grpSLOPE`), an ADMM solver has been implemented for the Group SLOPE model (function `admmSolverGroupSLOPE `).

# grpSLOPE 0.2.1

* Fix for the installation error on r-oldrel-windows-ix86+x86_64 (R v3.2.5) in the CRAN package check results. This error was caused by the generic S3 method `sigma()` not being available from the `stats` package prior to R v3.3.0.

* Checks for missing data in the inputs X, y and group were added in grpSLOPE().

* A check for whether the input matrix X has columns with 0 variance is performed in grpSLOPE() when normalize=TRUE. This should prevent division by 0, when the columns of X are standardized to have norms equal to 1.
