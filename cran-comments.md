## Test environments

* local Ubuntu and Windows, R (release)
* ubuntu 12.04 (on travis-ci), R (oldrel, release, devel)
* Windows (on AppVeyor), R (release)
* win-builder, R (release, devel)
* r-hub for UBSAN checks (Debian Linux, R-devel, GCC ASAN/UBSAN)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

There are no reverse dependencies.

## Comments

* Fixed memory access errors suggested bu UBSAN checks on CRAN
* Fixed spelling in the documentstion
