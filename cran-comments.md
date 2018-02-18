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

This is a minor update to make the package compatibile with the 
enchancements in the future package that it depends on.

This release fixes minor bug in the tests (#2) that can return
errors when testing the package on single core machine. I switched
to using future.apply package, since the future_lapply function
was moved there from the future package.

