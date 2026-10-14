# Validation of the documented development version

Development version: 1.0.0.9000. Local validation performed on 9 October 2026
using Windows and R 4.4.1, before this repository update.

* The regression suite passed: 19 cases and 82 assertions, including real
  classification engines and the installed optional FP-growth backend.
* `R CMD check --no-manual` finished with `Status: OK` for the corrected package
  and its then-current vignettes.
* After adding the complete English tutorial, the package build, installation
  and line-by-line English example succeeded. All 12 functions and all 128
  formal parameters are covered in that manual. The full package check was not
  repeated after these documentation additions.
* Both the reduced MSc example and synthetic raw-count example completed.
* Comparison with the archived MSc/article supplements is described in
  [RESULTS_COMPARISON_EN.txt](RESULTS_COMPARISON_EN.txt). It verifies specified
  outputs; it is not a complete rerun of the scientific study.

Runtime warnings included dependencies compiled under R 4.4.3, automatic numeric
discretization in rule mining and the explicit mining rule-length limit. The
manual explains how discretization and mining parameters affect interpretation.

The GitHub Actions workflow checks future pushes independently. Its result must
be read from the Actions tab; a local successful check does not imply that a
remote run has already passed.

## Run the examples

```r
library(MettleRNASeq)
source(system.file("examples", "run_example.R", package = "MettleRNASeq"))
run_mettle_example("mestrado", output_dir = "example_mestrado")
run_mettle_example("synthetic", output_dir = "example_synthetic")
```

For individual function calls in order, open
[`step_by_step.R`](../inst/examples/step_by_step.R) and use the
[complete English manual](MettleRNASeq_manual_completo_EN.html).
