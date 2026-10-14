<img src="man/figures/MettleRNASeq.png" alt="MettleRNASeq" height="100"/>

# MettleRNASeq

Select genes and explore relationships in RNA-seq expression data using
normalization, PCA, classification, correlations and association rules.

**Author and maintainer:** Samella Salles.

## Install

Requires R >= 4.2.0.

```r
install.packages(c("BiocManager", "remotes", "randomForest"))
BiocManager::install(c("edgeR", "DESeq2"), ask = FALSE, update = FALSE)
remotes::install_github("SamellaSalles/MettleRNASeq", build_vignettes = FALSE)
library(MettleRNASeq)
```

## Run an example

```r
source(system.file("examples", "run_example.R", package = "MettleRNASeq"))
result <- run_mettle_example("mestrado", output_dir = "example_results")
```

This example uses eight bone samples and 235 genes. The data are already
normalized. It saves tables, results, plots and session information in
`example_results`. For a synthetic example starting from raw counts:

```r
raw_result <- run_mettle_example("synthetic", output_dir = "synthetic_results")
```

## Use your own data

Provide expression values with **samples in rows and genes in columns**,
unique sample and gene names, and sample conditions in matching row order.
Normalize raw counts once; skip normalization for processed expression data.

Follow the [step-by-step manual](docs/MettleRNASeq_step_by_step_manual.html)
and [script](inst/examples/step_by_step.R) to run each function in order.
Download the manual using **Download raw file**, then open it in your browser.
The manual explains inputs, parameters, outputs and interpretation.

For individual function help:

```r
?normalize_count_matrix
?classification_approach_advanced
```

The examples use reduced settings. Classification scores are exploratory;
they do not establish performance on unseen patients. Rule graphs show gene
associations. The manual explains replicate partitions and discretization.

## Changes and support

See [NEWS.md](NEWS.md) for fixes and behavior changes. Checked MSc outputs
matched the archived values within numerical precision; the complete analysis
was not rerun. `icc_difference(..., unique_pairs = FALSE)` restores the
historical full-table layout. This is development version 1.0.0.9000.

[Report an issue](https://github.com/SamellaSalles/MettleRNASeq/issues) with a
small example and `sessionInfo()`. See [CONTRIBUTING.md](CONTRIBUTING.md) to
contribute.

## Citation and licence

Use `citation("MettleRNASeq")` for the software citation. For the method,
see the [preprint](https://doi.org/10.1101/2025.05.06.652387).
Example data attribution is in [PROVENANCE.txt](inst/extdata/PROVENANCE.txt).
Code is licensed under [GPL-2.0-or-later](COPYING).

## Acknowledgements

Otávio Brustolini, Luciane Ciapina and Kary Ocaña supervised the master's
research underlying the package. LNCC, CAPES and Rebecca Salles are acknowledged
for historical support and discussions.
