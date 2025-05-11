<img src="MettleRNASeq.png" alt="MettleRNASeq" class="inline" height="100"/>

# MettleRNASeq Package for R: Complex RNA-Seq Data Analysis and Gene Relationships Exploration Based on Machine Learning
__Current Version:__ 1.0
__Date:__ 2025-05
__Authors:__ Samella Salles (<ssalles@posgrad.lncc.br>), Otávio Brustolini (<tavinbio@lncc.br>), Luciane Ciapina (<luciane@lncc.br>), Kary Ocaña (<karyann@lncc.br>)

**Description**: 
#### This package follows the MettleRNASeq approach described in the bioRxiv preprint: https://doi.org/10.1101/2025.05.06.652387


**MettleRNASeq** is a robust alternative for RNA-Seq data analysis, particularly tailored for complex datasets that are common in biological research. Traditional differential gene expression (DGE) analysis might struggle when RNA-Seq datasets possess characteristics that hinder the power of statistical analyses, such as a limited number of replicates or high variability. MettleRNASeq integrates machine learning techniques, a tailored classification approach, association rule mining, and complementary correlation analysis to accurately identify key genes that distinguish experimental conditions and emphasize gene relationships.

This package is designed to be easy to use and practical, especially for biologists. It fits a wide range of transcriptomic studies, providing users with flexible tools to analyze complex datasets and gain insights into disease mechanisms, treatments, and their progression.
The package offers a range of functions for RNA-Seq data analysis, including:
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

__Acknowledgements:__ The authors thank the National Laboratory for Scientific Computing (LNCC, Brazil) for providing support and access to the Santos Dumont (SDumont) Supercomputer. We also acknowledge the Coordination for the Improvement of Higher Education Personnel (CAPES) for funding this research. Special thanks to Rebecca Salles for the valuable discussions and feedback.

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

#### **normalize_count_matrix** – Data Cleaning and Normalization

Normalize RNA-Seq count data using TMM, DESeq2, or scaling methods.

```r
# Example count matrix (samples as rows, genes as columns)
count_matrix <- data.frame(
  Gene1 = c(10, 20, 30, 40),
  Gene2 = c(15, 25, 35, 45),
  Gene3 = c(12, 22, 32, 42)
)

# Normalize using TMM method
norm_counts_tmm <- normalize_count_matrix(count_matrix, method = "TMM")

# Normalize using DESeq2 method
norm_counts_deseq2 <- normalize_count_matrix(count_matrix, method = "DESeq2")

# Normalize using scaling method
norm_counts_scale <- normalize_count_matrix(count_matrix, method = "scale")
```

#### **compute_distance_matrix** - Generation of distance matrices and heatmaps
Compute distance matrices and generate heatmaps for visualizing sample similarities.

```r
# Compute distance matrix using Euclidean distance and plot heatmap
dist_matrix_euclidean <- compute_distance_matrix(count_matrix, use_correlation = FALSE, method = "euclidean")

# Compute distance matrix using Spearman correlation and plot heatmap
dist_matrix_spearman <- compute_distance_matrix(count_matrix, use_correlation = TRUE, correlation_method = "spearman")

# Compute distance matrix without plotting heatmap
dist_matrix_no_plot <- compute_distance_matrix(count_matrix, plot_heatmap = FALSE)
```

#### **pca_analysis** – PCA Analysis and Plotting
Perform PCA analysis, visualize principal components, and generate 2D and 3D PCA plots.

```r
# Example 1: PCA analysis with a design matrix
  # Example count matrix (samples as rows, genes as columns)
 count_matrix <- data.frame(
   Gene1 = c(10, 20, 30, 40),
   Gene2 = c(15, 25, 35, 45),
   Gene3 = c(12, 22, 32, 42)
 )

# Example design matrix (conditions for each sample)
 design_matrix <- data.frame(
  Condition = c("A", "A", "B", "B")
 )

# Perform PCA analysis and generate plots
pca_results <- pca_analysis(count_matrix, design_matrix = design_matrix, color_by = "Condition", shape_by = NULL, center = TRUE, scale = TRUE, plot_3d = TRUE)
 
# Access PCA results
 pca_results$pca_result
 pca_results$eigenvalues
 pca_results$cos2

# Example 2: PCA analysis with color and shape vectors (no design matrix)
 color <- c("F", "F", "S", "S")
 shape <- c("FA", "FB", "SA","SB")

# Perform PCA analysis with custom color and shape vectors, generating 2D and 3D plots
 pca_results <- pca_analysis(count_matrix, plot_2d = TRUE, plot_cos2 = FALSE, design_matrix = NULL, color_by = color, shape_by = shape, plot_3d = TRUE)

# Perform PCA analysis and generate 2D and 3D plots
pca_results <- pca_analysis(count_matrix, plot_2d = TRUE, plot_3d = TRUE)

# Access PCA components
pca_results$pca_result  # PCA results object
pca_results$variance    # Variance explained by each component
pca_results$cos2        # Cos2 values for assessing representation quality
```

#### **feature_selection_pca** – Feature Selection Based on PCA
Select features based on PCA results, allowing users to choose the best and worst genes based on their contribution to variation from both the first and last principal components.

```r
# Perform PCA and feature selection based on the top genes of the first principal component
pca_result <- prcomp(count_matrix, center = TRUE, scale = TRUE)
feature_selection_results <- feature_selection_pca(pca_result, count_matrix, num_pcs = 2, num_genes_per_pc = 1)

# Access selected genes
feature_selection_results$selected_genes  # List of selected genes
feature_selection_results$selected_count_matrix  # New count matrix with selected genes
```

#### **prepare_input_matrix** – Preparing Input Matrix for Classification
Prepare input matrices for classification tasks by organizing data into a format suitable for machine learning algorithms.

```r
# Example sample and condition vectors
sample <- c("F1B", "F2B", "F3B", "F4B", "S1B", "S2B", "S3B", "S4B")
condition <- c("FLASH", "FLASH", "FLASH", "FLASH", "STAN", "STAN", "STAN", "STAN")

# Example count matrix
data <- matrix(runif(8 * 100), nrow = 100, ncol = 8)
colnames(data) <- sample

# Prepare input matrix for classification
result <- prepare_input_matrix(count_matrix = data, sample_vector = sample, condition_vector = condition)
result
```

#### **classification_approach_advanced** – Advanced Classification with Gene Resampling
Apply advanced classification techniques, including gene resampling and partitioning, to improve model accuracy and robustness.

```r
# Example classification setup with resampling and hyperparameters
hyper_params <- list(
   rf = data.frame(mtry = c(1, 5, 10)))

# Perform advanced classification with gene resampling
classification_results <- classification_approach_advanced(
   input_data = result,
   algorithms = c("rf"),
   hyper_params = hyper_params,
   num_resamples = 5,
   num_genes = 50
 )

# Access evaluation metrics and models
classification_results$eval_table    # Evaluation table with accuracy metrics
classification_results$models        # Trained models per algorithm and resample
```

#### **genes_from_classification_results** – Selecting Genes Based on Classification Results
Select genes based on classification results with a specified accuracy cutoff.

```r
# Select genes based on classification results with a specified accuracy cutoff
genes_results <- genes_from_classification_results(
  eval_table = classification_results$eval_table,
  genes_table = classification_results$genes_table,
  count_matrix = data,
  cutoff = 0.6,
  metric_column = "summary_test_accuracy",
  selection_mode = "unique",
  plot_type = "both"
)

# Access selected genes and filtered count matrix
genes_results$selected_gene_positions  # Selected gene positions (column indices)
genes_results$filtered_matrix          # Filtered count matrix with selected genes
```

#### **icc_analysis** – Intra-Condition Correlation Analysis
Perform intra-condition correlation analysis for gene relationships.

```r
# Perform intra-condition correlation analysis for gene relationships
icc_results <- icc_analysis(data, condition_split = c(2, 4), method = "spearman", plot = TRUE)

# Access correlation matrices for condition 1 and condition 2
icc_results$correlation_matrix_cond1  # Correlation matrix for condition 1
icc_results$correlation_matrix_cond2  # Correlation matrix for condition 2
```

#### **icc_difference** – Absolute Difference Between Intra-Condition Correlations
Compute the absolute difference between intra-condition correlations (ICC) to identify unique and shared gene relationships.

```r
# Compute absolute difference between intra-condition correlations
icc_diff_results <- icc_difference(icc_results, plot = TRUE)

# Access the results of the absolute difference
icc_diff_results$diff_matrix  # Matrix of absolute differences between conditions
```

#### **mine_gene_associations** – Gene Association Mining
Identify conditional gene relationships using association rule mining.

```r
# Perform association rule mining for gene relationships with different plot types and methods
association_results <- mine_gene_associations(data, 1:4, plot_type = "scatter", interactive = TRUE, top_n = 600)
association_results <- mine_gene_associations(data, 1:4, method = "fpgrowth", plot_type = "graph", interactive = FALSE, lift = 2)

# Focusing on genes of interest
association_results <- mine_gene_associations(data, 1:4, plot_type = "graph", interactive = TRUE, top_n = 10, target_genes = "Gene_1", search_type = "lhs")

# Access the identified gene associations
association_results$rules  # Association rules identified between genes
```

#### **compare_gene_rules** – Comparing Gene Association Rules
Compare gene association rules to highlight gene hubs that are unique or shared between experimental conditions.

```r
# Compare gene association rules between conditions
association_results <- mine_gene_associations(data, 1:4)
association_results_x <- mine_gene_associations(data, 5:8)
comparison_results <- compare_gene_rules(rules_df1 = association_results$rules_df, rules_df2 = association_results_x$rules_df,
                      condition_names = c("F", "S"), plot = TRUE)

# Access the comparison of gene hubs between conditions
comparison_results$gene_hubs  # Gene hubs unique or shared across conditions
```



