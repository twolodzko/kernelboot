## Test environments

* local Ubuntu and Windows, R (release)
* ubuntu 12.04 (on travis-ci), R (oldrel, release, devel)
* Windows (on AppVeyor), R (release)
* win-builder, R (release, devel)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

This is a new release, so there are no reverse dependencies.

## Comments

A number of tests is included in tests/ to ensure that things work as expected.
Moreover, a number of tests was conducted independently. Additional tests and
checks are done automatically outside of CRAN environment (are ignored for CRAN).
The Rcpp functions were benchmarked against the pure-R alternatives.
Multiple examples are provided to "visually" check the output. Extensive
documentation accompanies the package.
