# MettleRNASeq 1.1.0

* Retain default gene resampling, replicate-column partitions,
  internal training procedure, and gene identity predictor by default.
* Fix classification with one resample and varying replicate counts. Evaluate
  each held-out replicate independently using the repeated-expression
  test representation; report its accuracy and the partition mean.
* Match all bestTune parameters, record fit failures, and keep incomplete
  resample summaries missing. Stop on failed fits by default; an explicit
  continuation option returns diagnostic tables.
* Match genes across conditions and downstream matrices by gene ID. Preserve
  original sample IDs in a replicate mapping table.
* Keep matrix dimensions and surviving sample names during normalization.
  Preserve the original non-centered scaling convention and document it.
* Compare complete gene sets or expression-state itemsets independently of
  rule quality metrics. Preserve original rules and project every member into
  the gene graphs; handle empty graphs and large cluster counts.
* Plot the recorded classification summaries and separate algorithms. Align
  both axes of correlation matrices before computing differences; return unique
  undirected pairs by default, with a legacy full-table option.
* Validate input labels and dimensions, fix target-gene matching in rule mining,
  remove automatic dependency installation during analysis, and repair the
  documentation example. Add regression tests using real classification engines.
* Pass support and confidence correctly to the FP-growth wrapper and require
  its optional backend explicitly. Fix PCA plot metadata and Plotly namespaces.
* Add a runnable bone example with processed expression values, provenance,
  a synthetic raw-count example, a user vignette and contributor documentation.
* Load mining dependencies explicitly without attaching packages to the search
  path. Add a GitHub Actions check workflow for future repository pushes.
