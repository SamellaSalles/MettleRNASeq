<img src="man/figures/MettleRNASeq.png" alt="MettleRNASeq" height="100"/>

# MettleRNASeq

An R package for exploratory gene selection and relationship analysis in RNA-seq
expression matrices. The workflow combines normalization, PCA, gene resampling
with replicate partitions, correlations and association rule mining.

**Author and maintainer:** Samella Salles.
**Development version documented here:** 1.0.0.9000 (unreleased).
**Original method:** [MettleRNASeq preprint](https://doi.org/10.1101/2025.05.06.652387).

The package implements the author's MSc workflow. It uses established tools
including edgeR, DESeq2, caret and arules and returns intermediate tables and
identifiers. Its exploratory scores and gene relationships complement
statistical analysis; they are not validated predictions for unseen patients
or evidence of causality.

## Installation

```r
install.packages(c("BiocManager", "remotes", "randomForest"))
BiocManager::install(c("edgeR", "DESeq2"), ask = FALSE, update = FALSE)
remotes::install_github("SamellaSalles/MettleRNASeq", build_vignettes = TRUE)
library(MettleRNASeq)
packageVersion("MettleRNASeq")
```

Building the installed tutorials requires Pandoc, which is included with RStudio.
To install without building tutorials, use `build_vignettes = FALSE` and read
the documentation linked below.
R >= 4.2.0 is declared; this version was tested on Windows with R 4.4.1.
Other caret engines need their corresponding packages. FP-growth additionally
requires `fim4r`. Analysis functions do not install dependencies automatically.

## Reproducible example

The included MSc bone subset contains eight samples and 235 genes from GEO
GSE173944. It retains the author's original processed expression values.
It is **already normalized and preselected**, not raw counts. Its provenance
is documented in [PROVENANCE.txt](inst/extdata/PROVENANCE.txt).

```r
source(system.file("examples", "run_example.R", package = "MettleRNASeq"))
result <- run_mettle_example("mestrado", output_dir = "example_mestrado")
```

The script runs the full downstream workflow with three gene resamples,
exports CSV tables and RDS results, saves session information, and produces
a PDF of plots. A synthetic raw-count example also runs TMM normalization:

```r
raw_example <- run_mettle_example("synthetic", output_dir = "example_synthetic")
```

These short examples teach the software. They do not reproduce every analysis
in the dissertation or establish new biological results. Fixed seeds improve
repeatability within a recorded computational environment.

## Manuals and scientific interpretation

- [Complete English manual with figures and results](docs/MettleRNASeq_step_by_step_manual.html): download the HTML file using GitHub's **Download raw file** button and open it in a browser.
- [English step-by-step manual source](vignettes/MettleRNASeq-step-by-step.Rmd): every function, all parameter choices, actual results and links between stages.
- [English line-by-line script](inst/examples/step_by_step.R).
- [Runnable workflow](inst/examples/run_example.R).
- Installed function help: `help("classification_approach_advanced", package = "MettleRNASeq")`.
- Installed tutorial: `browseVignettes("MettleRNASeq")` after building vignettes.
- [Maintenance changes](NEWS.md).

Most functions expect samples in rows and genes in columns. Classification
uses a gene-condition representation with replicate columns, retains gene
identity by default, and preserves the original internal training procedure.
Read the manual before interpreting train or held-out replicate scores.
Automatic discretization is performed separately within each condition;
use explicitly shared state definitions when comparing expression-state items.
Rule graphs are gene projections, not direct regulatory networks.

## Validation against the MSc supplements

The checked bone and skin outputs matched the archived results within numerical
precision: normalization, distances, PCA summaries, PCA-selected genes, archived
gene resamples, correlations and comparison of archived mined rules. Eight
reconstructed random-forest partition fits also matched the original script
and archived scores. The full 400,000-fit search and full rule mining were not
rerun. The quick examples use smaller settings than the research analysis.

The ICC-difference table now lists each undirected off-diagonal pair once.
Use `icc_difference(..., unique_pairs = FALSE)` for the historical full-table
representation. Pairwise coefficients in the checked supplements are unchanged.
The full model search selecting the historical resamples, differential-expression
and enrichment analyses were not rerun. These checks do not constitute independent
biological validation.

The local regression suite passed 19 cases and 82 assertions on Windows with
R 4.4.1. A package check completed with `Status: OK`; after adding the complete
English tutorial, building, installing and executing that tutorial also succeeded.
Current remote checks are available in the repository's Actions tab.

![Original workflow](man/figures/Fluxogram_MettleRNASeq.png)

## Development and support

[Report an issue](https://github.com/SamellaSalles/MettleRNASeq/issues) with a
minimal example, version and `sessionInfo()`. See [CONTRIBUTING.md](CONTRIBUTING.md)
for development and support expectations. The repository includes regression
tests and a GitHub Actions package-check workflow; CI runs after a push.

Code is distributed under **GPL-2.0-or-later**. See [COPYING](COPYING).
The software citation metadata in `CITATION.cff` refers to the existing 1.0.0
release; update it when a new release is actually made.

## Acknowledgements

Otávio Brustolini, Luciane Ciapina and Kary Ocaña supervised the master's
research underlying the package. Historical support from LNCC and CAPES,
and discussions with Rebecca Salles, are acknowledged in the original project.
