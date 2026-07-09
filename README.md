<img src="man/figures/MettleRNASeq.png" alt="MettleRNASeq" class="inline" height="100"/>
# MettleRNASeq Package for R: Machine-learning-assisted analysis of complex RNA-seq datasets
__Current Version:__ 1.0.0
__Author:__ Samella Salles (<ssalles@posgrad.lncc.br>)
__Advisors:__ Otávio Brustolini (<tavinbio@lncc.br>), Luciane Ciapina (<luciane@lncc.br>), Kary Ocaña (<karyann@lncc.br>)

**Description**: 

**MettleRNASeq** is a robust alternative for RNA-Seq data analysis, particularly tailored for complex datasets that are common in biological research. Traditional differential gene expression (DGE) analysis might struggle when RNA-Seq datasets possess characteristics that hinder the power of statistical analyses, such as a limited number of replicates or high variability. MettleRNASeq integrates machine learning techniques, a tailored classification approach, association rule mining, and complementary correlation analysis to accurately identify key genes that distinguish experimental conditions and emphasize gene relationships.

This package is designed to be easy to use and practical, especially for biologists. It fits a wide range of transcriptomic studies, providing users with flexible tools to analyze complex datasets and gain insights into disease mechanisms, treatments, and their progression.

MettleRNASeq is accompanied by a methodological preprint (*described in the bioRxiv: https://doi.org/10.1101/2025.05.06.652387*) and provides configurable outputs for exploratory analysis, including classification metrics, resampling summaries, selected gene sets, correlation matrices, association rules, rule comparisons, and multiple visualization options. Several functions include user-defined parameters and optional parallel execution to support flexible transcriptomic analyses.

**Overview of the MettleRNASeq approach**
![MettleRNASeq workflow](man/figures/Fluxogram_MettleRNASeq.png)
#####Data obtention steps are shown in the white rectangles. MettleRNASeq approach is illustrated in the light gray box, encompassing data preparation steps (in light pink) and the three key analyses: tailored classification (in green) and the gene relationship-focused correlation analysis and association rule mining (in blue tones). Biological analyses both gene and relationship-focused are shown in yellow.

MettleRNASeq can be gene and relationship focused and presents as main features:

✔ Robust preprocessing

✔ PCA-based feature selection

✔ Machine-learning-assisted classification

✔ Correlation analysis

✔ Association Rule Mining

✔ Network visualization

✔ Candidate biomarker discovery


This package offers a range of functions for RNA-Seq data analysis, including:
*  **normalize_count_matrix** - Data cleaning and normalization using TMM, DESeq2, or scaling methods.
*  **compute_distance_matrix** - Generation of distance matrices and heatmaps.
*  **pca_analysis** - PCA analysis, visualization of principal components, plotting of cos2, 2D and 3D PCA plots.
*  **feature_selection_pca** - Feature selection based on PCA results, allowing choice of best and worst genes for variation explanation, from first and last PCs. 
*  **prepare_input_matrix** - Preparation of input matrices for classification.
*  **classification_approach_advanced** - Advanced classification techniques, including gene resampling and partitioning.
*  **genes_from_classification_results** - Filtering classification results for gene selection, plotting of classification results and updating the count matrix.
*  **acc_analysis** - Analysis of gene relationships through correlation-based approaches: Across-Condition Correlations
*  **icc_analysis** - Analysis of gene relationships through correlation-based approaches: Intra-Condition Correlations
*  **icc_difference** - Absolute Difference Between Intra-Condition Correlations (ICC) for identification of unique and shared gene relationships between conditions.
*  **mine_gene_associations** - Analysis of conditional gene relationships through association rule mining.
*  **compare_gene_rules** - Comparison of Gene Association Rules, identifying gene hubs unique and shared between conditions.

---
### Usage:
```r
#Install MettleRNASeq package
install.packages("devtools")
devtools::install_github("SamellaSalles/MettleRNASeq")

#Load MettleRNASeq package
library("MettleRNASeq")
```
#####


## Quick workflow example

The example below illustrates the main MettleRNASeq workflow using a small simulated RNA-seq count matrix with two experimental conditions and four replicates per condition.

```r
library(MettleRNASeq)

set.seed(123)

# ------------------------------------------------------------
# 1. Create an example RNA-seq count matrix
# ------------------------------------------------------------
# Rows = samples
# Columns = genes

sample_names <- c(
  paste0("FLASH_", 1:4),
  paste0("STANDARD_", 1:4)
)

condition <- c(
  rep("FLASH", 4),
  rep("STANDARD", 4)
)

count_matrix <- as.data.frame(matrix(
  rnbinom(8 * 80, size = 10, mu = 100),
  nrow = 8,
  ncol = 80
))

rownames(count_matrix) <- sample_names
colnames(count_matrix) <- paste0("Gene_", 1:80)

# Add a small condition-related signal to a subset of genes
count_matrix[condition == "FLASH", 1:10] <-
  count_matrix[condition == "FLASH", 1:10] + 50

# ------------------------------------------------------------
# 2. Normalize and inspect sample similarity
# ------------------------------------------------------------

norm_counts <- normalize_count_matrix(
  count_matrix = count_matrix,
  method = "TMM",
  min_count = 10,
  group_size = 4
)

distance_matrix <- compute_distance_matrix(
  count_matrix = norm_counts,
  use_correlation = FALSE,
  method = "euclidean",
  plot_heatmap = TRUE
)

# ------------------------------------------------------------
# 3. Explore data structure with PCA
# ------------------------------------------------------------

design_matrix <- data.frame(
  Sample = sample_names,
  Condition = condition
)

pca_results <- pca_analysis(
  count_matrix = norm_counts,
  design_matrix = design_matrix,
  color_by = "Condition",
  center = TRUE,
  scale = FALSE,
  plot_2d = TRUE,
  plot_3d = FALSE
)

# ------------------------------------------------------------
# 4. Select genes contributing most to PCA variation
# ------------------------------------------------------------

selected <- feature_selection_pca(
  pca_result = pca_results$pca_result,
  count_matrix = norm_counts,
  num_pcs = 2,
  num_genes_per_pc = 15
)

selected_genes <- selected$selected_genes
refined_matrix <- selected$selected_count_matrix

# ------------------------------------------------------------
# 5. Prepare the classification input
# ------------------------------------------------------------

classification_input <- prepare_input_matrix(
  count_matrix = refined_matrix,
  sample_vector = sample_names,
  condition_vector = condition
)

# ------------------------------------------------------------
# 6. Run the tailored classification approach
# ------------------------------------------------------------

hyper_params <- list(
  rf = data.frame(mtry = c(1, 3, 5))
)

classification_results <- classification_approach_advanced(
  input_data = classification_input,
  algorithms = c("rf"),
  hyper_params = hyper_params,
  num_resamples = 10,
  num_genes = 10,
  train_percentage = 0.75,
  summary_stat = "median",
  replace_genes = TRUE,
  parallel = FALSE,
  seed = 123
)

head(classification_results$eval_table)

# ------------------------------------------------------------
# 7. Select genes from the best-performing resamples
# ------------------------------------------------------------

genes_results <- genes_from_classification_results(
  eval_table = classification_results$eval_table,
  genes_table = classification_results$genes_table,
  count_matrix = refined_matrix,
  partition_table = classification_results$partition_table,
  cutoff = 0.6,
  metric_column = "summary_test_accuracy",
  selection_mode = "unique",
  plot_type = "both",
  hline_yintercept = 0.6
)

final_matrix <- genes_results$filtered_matrix

# ------------------------------------------------------------
# 8. Explore gene relationships
# ------------------------------------------------------------

# Across-condition correlations
acc_results <- acc_analysis(
  data = final_matrix,
  condition_split = c(4, 8),
  method = "spearman",
  plot = TRUE,
  orientation = "horizontal"
)

# Intra-condition correlations
icc_results <- icc_analysis(
  data = final_matrix,
  condition_split = c(4, 8),
  method = "spearman",
  plot = TRUE
)

# Absolute differences between intra-condition correlations
icc_diff <- icc_difference(
  icc_result = icc_results,
  plot = TRUE
)

# Association rule mining by condition
rules_flash <- mine_gene_associations(
  data = final_matrix,
  condition_rows = 1:4,
  method = "apriori",
  plot_type = "graph",
  interactive = FALSE,
  lift = 1
)

rules_standard <- mine_gene_associations(
  data = final_matrix,
  condition_rows = 5:8,
  method = "apriori",
  plot_type = "graph",
  interactive = FALSE,
  lift = 1
)

# Compare gene association rules between conditions
rule_comparison <- compare_gene_rules(
  rules_df1 = rules_flash$rules_df,
  rules_df2 = rules_standard$rules_df,
  condition_names = c("FLASH", "STANDARD"),
  layout_type = "linear",
  plot = TRUE
)
rule_comparison$shared_rules
rule_comparison$unique_rules_condition1
rule_comparison$unique_rules_condition2
```

This workflow follows the main logic of MettleRNASeq: normalize and inspect the data, identify informative genes, use the tailored classification approach to select discriminative genes, and then explore gene relationships through correlation analysis and association rule mining.


__Acknowledgements:__ The authors thank the National Laboratory for Scientific Computing (LNCC, Brazil) for providing support and access to the Santos Dumont (SDumont) Supercomputer. We also acknowledge the Coordination for the Improvement of Higher Education Personnel (CAPES) for funding this research. Special thanks to Rebecca Salles for the valuable discussions and feedback.
