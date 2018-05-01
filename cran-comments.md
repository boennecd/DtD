## Test environments
* `rhub::check_for_cran()`
* `rhub::check(platform = "linux-x86_64-rocker-gcc-san")`
* local Windows 10 machine with R 3.5.0
* win-builder (devel and release)
* Ubuntu 14.04 (on travis-ci), R 3.5.0
* Ubuntu 17.10 with clang 4.0.0, devel, and ASAN/UBSAN settings
* Ubuntu 17.10 with gcc 7.2.0, devel, and valgrind-3.13.0

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Benjamin Christoffersen <boennecd@gmail.com>'

New submission

Possibly mis-spelled words in DESCRIPTION:
  Merton's (11:49)

I have added R-core as a copyright holder as I use the `nmmin` function from 
`R_ext/Applic.h`.

## Re-submission
This is a re-submission. In this version I have:

* Added a reference for Merton's distance to default model in the 'Description' 
  field.
