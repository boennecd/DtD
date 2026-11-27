## New version submission

CRAN archived DtD 0.2.2 on 2025-12-25 after check findings remained unresolved.
DtD 0.2.3 uses R's default C++ standard and removes the unused `LazyData` field.

## Test environments

The local check used `--as-cran`.

* Debian 13 (x86_64), R 4.6.2, GCC and G++ 14.2.0.

The hosted checks used `--no-manual --as-cran` with incoming feasibility checks disabled.

* Ubuntu 26.04.1 LTS (x86_64), R 4.5.3, R 4.6.2, and R-devel r90998,
  GCC and G++ 15.2.0.
* macOS 26.7.2 (arm64), R 4.6.2, Apple Clang 21.0.0.
* Windows Server 2025 (x86_64), R 4.6.2 (UCRT), GCC and G++ 14.3.0.

## R CMD check results

`R CMD check --as-cran DtD_0.2.3.tar.gz` reported 0 errors, 0 warnings,
and 1 note. The incoming feasibility check records the previous CRAN archive.

All five hosted checks reported 0 errors, 0 warnings, and 0 notes. An earlier
R-devel check flagged deprecated `.Names` attributes in test fixtures. The
current fixtures use named vectors and lists.

## Sanitizer and Valgrind checks

All four checks used R-devel r90991 on x86_64 with `--no-manual`. Each check
ran the examples and passed all 506 test assertions.

* Fedora Linux 44, GCC and G++ 16.2.1. The check used ASan and UBSan with leak
  detection disabled. It used `--as-cran --no-build-vignettes`. It reported
  0 errors, 0 warnings, and 1 note because Ghostscript was unavailable for
  the PDF size check.
* Fedora Linux 44, GCC and G++ 16.2.1, with Valgrind 3.27.1. The check used
  `--use-valgrind --no-stop-on-test-error`, with origin tracking and full
  leak checks. It reported 0 errors, 0 warnings, and 0 notes. Valgrind
  reported 0 memory errors and 0 bytes definitely, indirectly, or possibly lost.
* Ubuntu 26.04.1 LTS, Clang 23.1.6. One check enabled ASan and UBSan with
  leak detection disabled. The other enabled UBSan. Both used
  `--no-stop-on-test-error` and reported 0 errors, 0 warnings, and 0 notes.
