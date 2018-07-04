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

Explicitely declaring dependencies for the unit tests.

As discussed during the discussion on CRAN@r-project.org,
I attempted to fix the clang-UBSAN warning, but I wasn't
able to confirm this warning, as I wasn't able to reproduce
it. Nonetheless, the potential cause was fixed.

