## Test environments

* local Arch Linux, R version 3.6.3 (2020-02-29), x86_64-pc-linux-gnu (64-bit)
* Docker container `rocker/r-base`, R version 3.6.3 (2020-02-29), x86_64-pc-linux-gnu (64-bit)
* Docker container `rocker/r-devel`, R Under development (unstable) (2020-04-05 r78150), x86_64-pc-linux-gnu (64-bit)
* Ubuntu 16.04.6 LTS (on travis-ci), R 3.6.2 (2017-01-27)
* win-builder (devel, release, and oldrelease)

## R CMD check results

* There was 1 NOTE on local Arch Linux: "Compilation used the following non-portable flag(s): '-march=x86-64'". But I believe this to be a local issue that won't affect the install on other systems.
* There were no ERRORs, WARNINGs, or NOTEs on the other test environments listed above.

## CRAN-related changes

None.
