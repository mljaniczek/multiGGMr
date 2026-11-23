## Submission

This is a new submission.

The package depends on 'ggmtools' (functions shared with the related
package 'spikeyglass'), which is submitted to CRAN separately and must be
accepted first.

The C++ sources in src/ include the mGHS sampler of Claudio Busatto
(GPL >= 2), adapted for this package; he is listed as contributor and
copyright holder in Authors@R and the Copyright field, and every change is
recorded in inst/mghs/UPSTREAM.md.

## Test environments

* local macOS (aarch64), R 4.5
* (to run before submission) win-builder (devel, release), R-hub (Linux,
  Windows, macOS, sanitizers)

## R CMD check results

0 errors | 0 warnings | 1 note

* New submission.
