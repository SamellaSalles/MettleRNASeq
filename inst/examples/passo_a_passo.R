# MettleRNASeq: execute os blocos na ordem, ou selecione linhas no RStudio.
# Codigo extraido do percurso principal do manual em portugues.
# A pasta de saida padrao e temporaria; escolha uma permanente para seu projeto.

# Etapa: load
library(MettleRNASeq)
packageVersion("MettleRNASeq")
expressao <- read.csv(system.file("extdata", "mestrado_bone_expression.csv",
  package = "MettleRNASeq"), row.names = 1, check.names = FALSE)
metadados <- read.csv(system.file("extdata", "mestrado_bone_metadata.csv",
  package = "MettleRNASeq"), stringsAsFactors = FALSE)
stopifnot(identical(rownames(expressao), metadados$Sample))
dim(expressao)
metadados
expressao[1:3, 1:4]

# Etapa: normalize
set.seed(123)
contagens <- as.data.frame(matrix(rnbinom(8 * 60, size = 10, mu = 100), 8))
rownames(contagens) <- paste0("Amostra_", 1:8)
colnames(contagens) <- paste0("Gene_", 1:60)
normalizado_sintetico <- normalize_count_matrix(contagens, method = "TMM",
  min_count = 10, group_size = 4, remove_low_variance = TRUE)
dim(normalizado_sintetico)
attr(normalizado_sintetico, "removed_genes")
attr(normalizado_sintetico, "removed_samples")
# Continuamos o percurso real com os valores ja normalizados:
normalizada <- expressao

# Etapa: distance
distancias <- compute_distance_matrix(normalizada,
  use_correlation = TRUE, correlation_method = "spearman",
  plot_heatmap = TRUE, display_numbers = TRUE,
  main_title = "Distancia entre amostras: 1 - Spearman")
round(distancias[1:4, 1:4], 3)

# Etapa: pca
pca <- pca_analysis(normalizada, design_matrix = metadados,
  color_by = "Condition", center = TRUE, scale = FALSE,
  plot_2d = TRUE, plot_3d = FALSE)
data.frame(PC = seq_along(pca$variance),
  variancia_percentual = pca$variance,
  acumulada_percentual = pca$cumulative_variance)
head(pca$cos2)

# Etapa: features
selecao_pca <- feature_selection_pca(pca$pca_result, normalizada,
  num_pcs = 2, num_genes_per_pc = 12,
  select_worst = FALSE, use_last_pcs = FALSE)
length(selecao_pca$selected_genes)
head(selecao_pca$selected_genes)
refinada <- selecao_pca$selected_count_matrix
dim(refinada)

# Etapa: prepare
entrada <- prepare_input_matrix(refinada, design_table = metadados,
  sample_col = "Sample", condition_col = "Condition")
head(entrada)
attr(entrada, "sample_map")

# Etapa: classification
classificacao <- classification_approach_advanced(entrada,
  algorithms = "rf", hyper_params = list(rf = data.frame(mtry = c(1, 3))),
  num_resamples = 3, num_genes = 10, train_percentage = .75,
  metrics = "Accuracy", summary_stat = "median", replace_genes = TRUE,
  parallel = FALSE, seed = 123, include_gene = TRUE,
  error_handling = "stop")
classificacao$partition_table
head(classificacao$eval_table)
classificacao$status_table
head(classificacao$gene_ids_table)

# Etapa: select
escores <- unique(classificacao$eval_table[,
  c("algorithm", "resample", "summary_test_accuracy")])
escores
id_escolhido <- escores$resample[which.max(escores$summary_test_accuracy)]
selecao <- genes_from_classification_results(classificacao$eval_table,
  classificacao$genes_table, refinada, resamples = id_escolhido,
  gene_ids_table = classificacao$gene_ids_table,
  partition_table = classificacao$partition_table, selection_mode = "unique",
  plot_type = "both", plot_resamples_only = FALSE, return_plot = TRUE)
selecao$selected_resamples
selecao$selected_gene_ids
selecionada <- selecao$filtered_matrix
dim(selecionada)

# Etapa: acc
stopifnot(identical(metadados$Condition,
  rep(c("FLASH", "STANDARD"), each = 4)))
divisao <- c(4, 8)
acc <- acc_analysis(selecionada, condition_split = divisao,
  method = "spearman", plot = TRUE, orientation = "horizontal")
acc$results

# Etapa: icc
icc <- icc_analysis(selecionada, condition_split = divisao,
  method = "spearman", plot = TRUE, plot_type = "lower",
  show_coef = TRUE, tl.cex = .6, number.cex = .6, order = "alphabet")
round(icc$correlation_matrix_cond1[1:4, 1:4], 2)
round(icc$correlation_matrix_cond2[1:4, 1:4], 2)

# Etapa: differences
diferencas <- icc_difference(icc, cutoff = NULL,
  unique_pairs = TRUE, plot = TRUE, fontsize = 7)
head(diferencas$difference_table)
dim(diferencas$absolute_diff_matrix)

# Etapa: difference-filter
mudancas_grandes <- icc_difference(icc, cutoff = 1, plot = FALSE)
head(mudancas_grandes$difference_table)

# Etapa: rules
regras_flash <- mine_gene_associations(selecionada, condition_rows = 1:4,
  method = "apriori", supp = .5, conf = .8, minlen = 2, maxlen = 2,
  plot_type = "scatter", interactive = FALSE, top_n = 20)
regras_standard <- mine_gene_associations(selecionada, condition_rows = 5:8,
  method = "apriori", supp = .5, conf = .8, minlen = 2, maxlen = 2,
  plot_type = "scatter", interactive = FALSE, top_n = 20)
length(regras_flash$rules)
length(regras_standard$rules)
head(regras_flash$rules_df)
head(regras_standard$rules_df)
head(arules::itemLabels(regras_flash$transactions))

# Etapa: compare
comparacao <- compare_gene_rules(regras_flash$rules_df,
  regras_standard$rules_df, condition_names = c("FLASH", "STANDARD"),
  comparison = "genes", layout_type = "linear", plot = TRUE)
comparacao$shared_rules
head(comparacao$unique_rules_condition1)
head(comparacao$unique_rules_condition2)

# Etapa: compare-items
comparacao_itens <- compare_gene_rules(regras_flash$rules_df,
  regras_standard$rules_df, condition_names = c("FLASH", "STANDARD"),
  comparison = "items", plot = FALSE)
nrow(comparacao_itens$shared_rules)

# Etapa: save
pasta <- file.path(tempdir(), "MettleRNASeq_passo_a_passo")
dir.create(pasta, recursive = TRUE, showWarnings = FALSE)
saveRDS(list(normalizada = normalizada, metadados = metadados,
  distancias = distancias, pca = pca, selecao_pca = selecao_pca,
  entrada = entrada, classificacao = classificacao, selecao = selecao,
  acc = acc, icc = icc, diferencas = diferencas,
  regras_flash = regras_flash, regras_standard = regras_standard,
  comparacao = comparacao), file.path(pasta, "percurso_completo.rds"))
write.csv(classificacao$eval_table, file.path(pasta, "classificacao.csv"),
  row.names = FALSE)
write.csv(classificacao$gene_ids_table, file.path(pasta, "genes_sorteados.csv"),
  row.names = FALSE)
write.csv(selecionada, file.path(pasta, "expressao_selecionada.csv"))
write.csv(acc$results, file.path(pasta, "acc.csv"), row.names = FALSE)
write.csv(icc$correlation_matrix_cond1, file.path(pasta, "icc_FLASH.csv"))
write.csv(icc$correlation_matrix_cond2, file.path(pasta, "icc_STANDARD.csv"))
write.csv(diferencas$difference_table, file.path(pasta, "diferencas_icc.csv"),
  row.names = FALSE)
write.csv(regras_flash$rules_df, file.path(pasta, "regras_FLASH.csv"),
  row.names = FALSE)
write.csv(regras_standard$rules_df, file.path(pasta, "regras_STANDARD.csv"),
  row.names = FALSE)
write.csv(comparacao$shared_rules, file.path(pasta, "regras_compartilhadas.csv"),
  row.names = FALSE)
write.csv(comparacao$unique_rules_condition1,
  file.path(pasta, "regras_exclusivas_FLASH.csv"), row.names = FALSE)
write.csv(comparacao$unique_rules_condition2,
  file.path(pasta, "regras_exclusivas_STANDARD.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(pasta, "sessionInfo.txt"))
list.files(pasta)
