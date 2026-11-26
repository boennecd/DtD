## New version submission

CRAN archived DtD 0.2.2 on 2025-12-25 after check findings remained unresolved.
DtD 0.2.3 uses R's default C++ standard and removes the unused `LazyData` field.

## Test environments

The local check used `--as-cran`.

* Debian 13 (x86_64), R 4.6.2, GCC and G++ 14.2.0.

The hosted checks used `--no-manual --as-cran` with incoming feasibility checks disabled.

* Ubuntu 26.04.1 LTS (x86_64), R 4.5.3, R 4.6.2, and R-devel r90991,
  GCC and G++ 15.2.0.
* macOS 26.7.2 (arm64), R 4.6.2, Apple Clang 21.0.0.
* Windows Server 2025 (x86_64), R 4.6.2 (UCRT), GCC and G++ 14.3.0.

## R CMD check results

`R CMD check --as-cran DtD_0.2.3.tar.gz` reported 0 errors, 0 warnings,
and 1 note. The incoming feasibility check records the previous CRAN archive.

The initial hosted checks on R 4.5.3 and R 4.6.2 reported no errors, warnings,
or notes. The initial R-devel check reported one note and six test warnings
about deprecated `.Names` attributes in the test fixtures. This version uses
named vectors and lists for those fixtures.
