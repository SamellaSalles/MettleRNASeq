#Installing packages:----------
library(plotly)
library(ggplot2)
library(DESeq2)
library(edgeR)
library(plyr)
library(tidyverse)
library(dplyr)
library(janitor)
library(tidyr)
library(corrplot)
library(pheatmap) 
library(RColorBrewer)
library(factoextra)
library(stringr)
library(caret)
library("rsample")
library(doParallel)  # Load doParallel package
library("dplyr")
library("tidyselect")
library(arules)
library(arulesViz)

#With the count matrix obtained after HTSeq Counts:
#Already in the .RDATA (count_matrix_B) -> Create a single dataframe with the count matrices!!!!! ------------
#------- #Criando uma lista dos meus dados para trocar nomes das colunas
#A coluna vai alterar de acordo com o nome do dataframe

# Function to rename columns in a data frame based on file name
rename_columns <- function(df, file_name) {
  # Generate new column names
  new_names <- c("Gene", paste0(file_name))
  # Rename columns
  names(df) <- new_names
  return(df)
}

# List of data frames
list_of_dfs <- list(F1B.merge, F2B.merge, F3B.merge, F4B.merge,
                    S1B.merge, S2B.merge, S3B.merge, S4B.merge,
                    U1B.merge, U2B.merge, U3B.merge, U4B.merge
                    #,F1S.merge, F2S.merge, F3S.merge, F4S.merge,
                    #S1S.merge, S2S.merge, S3S.merge, S4S.merge,
                    #U1S.merge, U2S.merge, U3S.merge, U4S.merge
                    )

# List of file names corresponding to each data frame
file_names <- c("F1B", "F2B", "F3B", "F4B",
                "S1B", "S2B", "S3B", "S4B",
                "U1B", "U2B", "U3B", "U4B"
                #,"F1S", "F2S", "F3S", "F4S",
               # "S1S", "S2S", "S3S", "S4S",
                #"U1S", "U2S", "U3S", "U4S"
               )

# Apply renaming function to each data frame in the list
for (i in seq_along(list_of_dfs)) {
  list_of_dfs[[i]] <- rename_columns(list_of_dfs[[i]], file_names[i])
}

# Accessing the first data frame to check the changes
print(head(list_of_dfs[[1]]))



----- #Juntar dataframes em 1 só
listdf <- list_of_dfs
listdf <- join_all(listdf, by = "Gene", type = "full", match = "all")
#There is no NA in the dataframe:
#sum(is.na(listdf))

-------- #Preparar data frame
count_matrix_B <- listdf
count_matrix_B <- t(listdf)
count_matrix_B <- count_matrix_B %>%
  row_to_names(1, remove_row = TRUE, remove_rows_above = TRUE) 
count_matrix_B <- count_matrix_B[, -c((ncol(count_matrix_B) - 4):ncol(count_matrix_B))]

# Convert "count_matrix_B" to a dataframe if it's not already
if (!is.data.frame(count_matrix_B)) {
  count_matrix_B <- as.data.frame(count_matrix_B)
}
# Convert character columns to numeric
for (col in names(count_matrix_B)) {
  count_matrix_B[, col] <- as.numeric(count_matrix_B[, col], na.rm = TRUE)
}
#count_matrix_B <-  sapply(count_matrix_B, as.numeric)
#return the datatype of each column
#head(print(sapply(count_matrix_B, class)))
#There is no NA in the dataframe:
#sum(is.na(count_matrix_B))

#Remover colunas sem vari?ncia nos valores entre os genes
#cols_to_keep <- sapply(count_matrix_B, function(x) var(x) != 0)
#count_matrix_B <- count_matrix_B[, cols_to_keep]
#PRONTO!

# ---------Normalize the data first (columns GENE): ----------------
  #DESeq2 with TMM normalization from edgeR

samples_B4 <- c("F1B", "F2B", "F3B", "F4B", "S1B", "S2B", "S3B", "S4B", "U1B", "U2B", "U3B", "U4B")
#samples_S4 <- c("F1S", "F2S", "F3S", "F4S", "S1S", "S2S", "S3S", "S4S", "U1S", "U2S", "U3S", "U4S")

# Define a function to extract condition and tissue
condition_map <- c(F = "FLASH", S = "Standard", U = "Untreated")

extract_condition <- function(sample_name) {
  condition <- condition_map[substr(sample_name, 1, 1)]
  return(paste(condition, sep = ""))
}
condition_B4 <- sapply(samples_B4, extract_condition)
#condition_S4 <- sapply(samples_S4, extract_condition)

directory <- "C:/Users/samel/OneDrive/Trabalho Pesquisa/Pesquisa Cancer HPC - Mestrado/Artigos/Paper Small Data Oxford BioInf/Code/input"
samplesFiles_B4 <- c("F1B.merge.fastq", "F2B.merge.fastq", "F3B.merge.fastq", "F4B.merge.fastq", "S1B.merge.fastq", "S2B.merge.fastq", "S3B.merge.fastq", "S4B.merge.fastq", "U1B.merge.fastq", "U2B.merge.fastq", "U3B.merge.fastq", "U4B.merge.fastq")
#samplesFiles_S4 <- c("F1S.merge.fastq", "F2S.merge.fastq", "F3S.merge.fastq", "F4S.merge.fastq", "S1S.merge.fastq", "S2S.merge.fastq", "S3S.merge.fastq", "S4S.merge.fastq", "U1S.merge.fastq", "U2S.merge.fastq", "U3S.merge.fastq", "U4S.merge.fastq")

sampleTable_B4 <- data.frame(row.names = NULL, sampleName = samples_B4, fileName = samplesFiles_B4, condition = condition_B4, stringsAsFactors = FALSE) # tissue = tissue
#sampleTable_S4 <- data.frame(row.names = NULL, sampleName = samples_S4, fileName = samplesFiles_S4, condition = condition_S4, stringsAsFactors = FALSE) # tissue = tissue

ddsHTSeq_B4 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable_B4,
                                          directory = directory,
                                          design= ~ condition) #group
ddsHTSeq_B4

#ddsHTSeq_S4 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable_S4,
 #                                         directory = directory,
  #                                        design= ~ condition) #group
#ddsHTSeq_S4

#Removing very low counts
smallestGroupSize <- 4 #Ideally the number of smallest samples you have for this condition/group
keep <- rowSums(counts(ddsHTSeq_B4) >= 10) >= smallestGroupSize
ddsHTSeq_B4 <- ddsHTSeq_B4[keep,]

ddsHTSeq_B4

# smallestGroupSize <- 4
# keep <- rowSums(counts(ddsHTSeq_S4) >= 10) >= smallestGroupSize
# ddsHTSeq_S4 <- ddsHTSeq_S4[keep,]
# 
# ddsHTSeq_S4

#B4 (with all samples from BONE)
tmm_data_B4 <- count_matrix_B %>% filter(!row_number() %in% 13:24) 
#Convert to numeric without losing rownames:
tmm_data_B4 <- data.frame(lapply(tmm_data_B4, function(x) as.numeric(as.character(x))),
                          check.names=F, row.names = rownames(tmm_data_B4))
tmm_data_B4 <- t(tmm_data_B4)
#colnames(tmm_SF_B_data)
#Rearranging the positions of the columns:
#tmm_SF_B_data <- tmm_SF_B_data[, c(5, 6, 7, 8, 1, 2, 3, 4)]

smallestGroupSize <- 4
keep <- rowSums(tmm_data_B4 >= 10) >= smallestGroupSize
tmm_data_B4 <- tmm_data_B4[keep,]

tmm_data_B4

tmm_B4 <- edgeR::calcNormFactors(tmm_data_B4, method="TMM")
# libsize is the column sum of the raw counts
N_B4 <- colSums(tmm_data_B4) #(vector of library sizes)
# size factor combines libsize with normalization factor into a single value
tmm.counts_B4 <- N_B4*tmm_B4/exp(mean(log(N_B4*tmm_B4)))

#Changing the sizefactor from DESeq2 (RLE) to TMM:
sizeFactors(ddsHTSeq_B4) <- tmm.counts_B4

ddsHTSeq_B4 <- DESeq(ddsHTSeq_B4)

res_B4_SxU <- results(ddsHTSeq_B4, contrast = c("condition", "Standard", "Untreated"), alpha = 0.05)#, alpha = 0.01
res_B4_FxU <- results(ddsHTSeq_B4, contrast = c("condition", "FLASH", "Untreated"), alpha = 0.05)#, alpha = 0.01
res_B4_SxF <- results(ddsHTSeq_B4, contrast = c("condition", "Standard", "FLASH"), alpha = 0.05)#, alpha = 0.01


#Getting the counts of the data that are already normalized (3 Conditions (3C)):
#Rows are GENES!
norm_B4_counts_3C <- counts(ddsHTSeq_B4, normalized=TRUE)

#Now rows are SAMPLES!
trans_norm_B4_counts_3C <- t(norm_B4_counts_3C)

#Same data but with just FLASH and STANDARD (2 Conditions (2C)):
norm_B4_counts <- norm_B4_counts_3C[1:8,]
trans_norm_B4_counts <- trans_norm_B4_counts_3C[1:8,]


#------------Similaridade entre amostras e condições--------------
#dist mostra distancia entre pares de LINHAS (rows)!
sim_B4 <- dist(trans_norm_B4_counts_3C,method="euclidean")
sim_B4 <- as.data.frame(as.matrix(sim_B4))
#COrrelação feita nas colunas, que aqui serão amostras!
#spearman_sim_B4 <- as.dist(1-cor(norm_B4_counts, method = "spearman"))
#spearman_sim_B4 <- as.data.frame(as.matrix(spearman_sim_B4))

num_breaks <- 5
colors <- colorRampPalette(c("#05287A", "#DCEEF0"))(num_breaks)

pheatmap(sim_B4,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         method = "complete",
         col=colors,
         breaks = seq(0, 500000, length.out = num_breaks))

#---------Correlação entre colunas (amostras) - Pode ser pesado computacionalmente!----------------
cor_B4 <- cor(norm_B4_counts_3C)

corrplot(cor_B4, method = "color", type = "lower", tl.cex = 0.5, tl.col="black",  number.cex = 0.5) #To add the numbers inside: addCoef.col = "black", #To add which values have significance: p.mat = cor_test_B4$p, sig.level = 0.05,


#-------------PCA to investigate the data (3 conditions (3C)):---------
pca_norm_B4_3C <- prcomp(trans_norm_B4_counts_3C, center = TRUE, scale = FALSE) # center = TRUE ,scale = FALSE
summary(pca_norm_B4_3C)
eig.val_norm_B4_3C <- get_eigenvalue(pca_norm_B4_3C)
eig.val_norm_B4_3C
#Investigate the cos2:
#Colored by cos2, meaning how much that sample is being well represented by the PCA
#Mais longe do centro, melhor representada é a amostra! Mais perto do centro, menos representada!

fviz_pca_ind(pca_norm_B4_3C,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#510000", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

#Correlation of how well the samples are being represented for each PCA
res.ind_norm_B4_3C <- get_pca_ind(pca_norm_B4_3C)
corrplot(res.ind_norm_B4_3C$cos2, is.corr=FALSE)

#PCA Plot
pca.data_norm_B4_3C <- data.frame(Sample=rownames(pca_norm_B4_3C$x),
                               X=pca_norm_B4_3C$x[,1],
                               Y=pca_norm_B4_3C$x[,2],
                               Condition=c("FLASH", "FLASH","FLASH","FLASH","STANDARD","STANDARD","STANDARD","STANDARD","UNTREATED","UNTREATED","UNTREATED","UNTREATED"))
pca.var_norm_B4_3C <- pca_norm_B4_3C$sdev^2
pca.var.per_norm_B4_3C <- round(pca.var_norm_B4_3C/sum(pca.var_norm_B4_3C)*100, 1)
ggplot(data=pca.data_norm_B4_3C, aes(x=X, y=Y, colour=Condition)) + #label=Sample,
  geom_point(size = 3) +
  xlab(paste("PC1 - ", pca.var.per_norm_B4_3C[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per_norm_B4_3C[2], "%", sep="")) +
  # scale_x_continuous(limits = c(4, 8), breaks = seq(4, 8, by = 1)) +
  #scale_y_continuous(limits = c(2, 7), breaks = seq(2, 7, by = 1)) +
  scale_color_manual(values = c("#0073e6", "#f57600", "#5ba300")) +
  scale_y_continuous(limits = c(-200000, 100000), breaks = seq(-200000, 100000, 100000)) +
  scale_x_continuous(limits = c(-200000, 300000), breaks = seq(-200000, 300000, 100000)) +
  theme_bw()


#----------#3D PCA Plots

pca_plotly_B4_3C <- plot_ly(as.data.frame(pca_norm_B4_3C$x), x = ~PC1, y = ~PC2, z = ~PC3, color = ~sampleTable_B4$condition) %>% add_markers()
pca_plotly_B4_3C

#------After elimination of UNTREATED condition, NEW PCAs are done for feature selection: ---------
pca_norm_B4 <- prcomp(trans_norm_B4_counts, center = TRUE, scale = FALSE) # center = TRUE ,scale = FALSE
summary(pca_norm_B4)
eig.val_norm_B4 <- get_eigenvalue(pca_norm_B4)
eig.val_norm_B4
#Investigate the cos2:
#Colored by cos2, meaning how much that sample is being well represented by the PCA
#Mais longe do centro, melhor representada é a amostra! Mais perto do centro, menos representada!

fviz_pca_ind(pca_norm_B4,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#510000", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

#Correlation of how well the samples are being represented for each PCA
res.ind_norm_B4 <- get_pca_ind(pca_norm_B4)
corrplot(res.ind_norm_B4$cos2, is.corr=FALSE)

#PCA Plot
pca.data_norm_B4 <- data.frame(Sample=rownames(pca_norm_B4$x),
                                  X=pca_norm_B4$x[,1],
                                  Y=pca_norm_B4$x[,2],
                                  Condition=c("FLASH", "FLASH","FLASH","FLASH","STANDARD","STANDARD","STANDARD","STANDARD"))
pca.var_norm_B4 <- pca_norm_B4$sdev^2
pca.var.per_norm_B4 <- round(pca.var_norm_B4/sum(pca.var_norm_B4)*100, 1)
ggplot(data=pca.data_norm_B4, aes(x=X, y=Y, colour=Condition)) + #label=Sample,
  geom_point(size = 3) +
  xlab(paste("PC1 - ", pca.var.per_norm_B4[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per_norm_B4[2], "%", sep="")) +
  scale_color_manual(values = c("#0073e6", "#f57600", "#5ba300")) +
  #scale_y_continuous(limits = c(-200000, 100000), breaks = seq(-200000, 100000, 100000)) +
  #scale_x_continuous(limits = c(-200000, 300000), breaks = seq(-200000, 300000, 100000)) +
  theme_bw()


#----------#3D PCA Plots
sampleTable_B4_2C <- sampleTable_B4[1:8,]
pca_plotly_B4 <- plot_ly(as.data.frame(pca_norm_B4$x), x = ~PC1, y = ~PC2, z = ~PC3, color = ~sampleTable_B4_2C$condition) %>% add_markers()
pca_plotly_B4


#--------------- Input refinement (feature selection): Top genes for each PC in PCA BONE 4 ------------------
## get the name of the top 100 measurements (genes) that contribute
## most to pc1.
loading_scores_B4 <- pca_norm_B4$rotation[,1]
gene_scores_B4 <- abs(loading_scores_B4) ## get the magnitudes
gene_score_ranked_PC1_B4 <- sort(gene_scores_B4, decreasing=TRUE)
top_100_genes_PC1_B4 <- names(gene_score_ranked_PC1_B4[1:100])

top_100_genes_PC1_B4 ## show the names of the top 10 genes

pca_norm_B4$rotation[top_100_genes_PC1_B4,1] ## show the scores (and +/- sign)

loading_scores_B4 <- pca_norm_B4$rotation[,2]
gene_scores_B4 <- abs(loading_scores_B4) ## get the magnitudes
gene_score_ranked_PC2_B4 <- sort(gene_scores_B4, decreasing=TRUE)
top_100_genes_PC2_B4 <- names(gene_score_ranked_PC2_B4[1:100])

top_100_genes_PC2_B4 ## show the names of the top 10 genes

pca_norm_B4$rotation[top_100_genes_PC2_B4,2] ## show the scores (and +/- sign)

loading_scores_B4 <- pca_norm_B4$rotation[,3]
gene_scores_B4 <- abs(loading_scores_B4) ## get the magnitudes
gene_score_ranked_PC3_B4 <- sort(gene_scores_B4, decreasing=TRUE)
top_100_genes_PC3_B4 <- names(gene_score_ranked_PC3_B4[1:100])

top_100_genes_PC3_B4 ## show the names of the top 10 genes

pca_norm_B4$rotation[top_100_genes_PC3_B4,3] ## show the scores (and +/- sign)

loading_scores_B4 <- pca_norm_B4$rotation[,4]
gene_scores_B4 <- abs(loading_scores_B4) ## get the magnitudes
gene_score_ranked_PC4_B4 <- sort(gene_scores_B4, decreasing=TRUE)
top_100_genes_PC4_B4 <- names(gene_score_ranked_PC4_B4[1:100])

top_100_genes_PC4_B4 ## show the names of the top 10 genes

pca_norm_B4$rotation[top_100_genes_PC4_B4,4] ## show the scores (and +/- sign)

#Obtendo os genes ?nicos de interesse dentro de lncRNA
top_genes_B4 <- c(top_100_genes_PC1_B4,top_100_genes_PC2_B4,top_100_genes_PC3_B4,top_100_genes_PC4_B4)
top_genes_B4 <- unique(top_genes_B4)

top_genes_B4 <- trans_norm_B4_counts[,top_genes_B4]

#--------Classification------------
ML_B4 <- as.data.frame(top_genes_B4) #or norm_B4_counts or raw_B4
ML_B4_F <- ML_B4 %>% filter(!row_number() %in% 5:8) #5:12
ML_B4_S <- ML_B4 %>% filter(!row_number() %in% 1:4) # %>% filter(!row_number() %in% 5:8)
#ML_B4_U <- ML_B4 %>% filter(!row_number() %in% 1:8)

ML_B4_F <- t(ML_B4_F) 
ML_B4_F <- as.data.frame(ML_B4_F)
names_F <- c( "B1" = "F1B", "B2" = "F2B", "B3" = "F3B", "B4" = "F4B")
ML_B4_F <- rename(ML_B4_F, all_of(names_F))

ML_B4_S <- t(ML_B4_S) 
ML_B4_S <- as.data.frame(ML_B4_S)
names_S <-  c( "B1" = "S1B", "B2" = "S2B", "B3" = "S3B", "B4" = "S4B")
ML_B4_S <- rename(ML_B4_S, all_of(names_S))

#ML_B4_U <- t(ML_B4_U) 
#ML_B4_U <- as.data.frame(ML_B4_U)
#names_U <- c( B1 = "U1B", B2 = "U2B", B3 = "U3B", B4 = "U4B")
#ML_B4_U <- rename(ML_B4_U, all_of(names_U))

ML_B4_F <- ML_B4_F %>%
  mutate(Condition = "FLASH")
ML_B4_S <- ML_B4_S %>%
  mutate(Condition = "STANDARD")
#ML_B4_U <- ML_B4_U %>%
#  mutate(Condition = "UNTREATED")
row_names_B4 <- as.data.frame(rownames(ML_B4_F) )
row_names_B4 <- bind_rows(row_names_B4,row_names_B4)#,row_names_B4
row_names_B4$Genes <- row_names_B4$`rownames(ML_B4_F)`
B4 <- bind_rows(ML_B4_F,ML_B4_S)#,ML_B4_U
B4 <- B4 %>%  mutate(Gene = row_names_B4$Genes)
rownames(B4) <- NULL
B4 <- B4[, c(6, 1, 2, 3, 4, 5)]

B4$Condition <- factor(B4$Condition)

#Using the gene refined after classification: In this case there were 235 genes!
B4_F <- B4[1:235,1:6] #<----- (top 100), B4[1:117,1:6] (top 50),  B4[1:57,1:6] (top 25), B4[1:24,1:6] (top 10),  B4[1:2236,1:6](top 1000), ou B4_top_genes_imp[1:55,1:6] ou B4_all_genes[1:21590,1:6] ou B4[1:2236,1:6] (for top 1000 genes)
data <- B4_F

#-----------------New way to get the model and genes important

calculate_accuracy <- function(data, indices, cols) {
  train_data <- data[indices, c(1,cols,6)]
  test_data <- data[indices,!(1:length(data)) %in% c(cols,6)]
  test_data <- cbind(test_data,replicate(length(cols)-1,test_data[2]))
  test_data$Condition <- data[indices,"Condition"]
  colnames(test_data) <- c("Gene",names(data[cols]),"Condition")
  indices_train_data <- seq_along(indices) #todos os indices em train_data (bootstrapinho)
  set.seed(123)
  model <- train(Condition~.,  
                 data=train_data,  
                 method="rf",
                 metric = "Accuracy", 
                 preProcess = "center",
                 trControl=trainControl(method="boot", number=1,index = list(indices_train_data),indexOut = list(indices_train_data)),
                 tuneLength = 10,
                 #tuneGrid = data.frame(mtry=mtry),
                 allowParallel= FALSE)
  acc_treino <- model$results%>%
    filter(mtry==model$bestTune$mtry)%>%
    select(Accuracy)
  predictions <- predict(model, newdata = test_data[1:5])
  accuracy <- mean(predictions == test_data$Condition)
  return(list(model=model,acc_treino=acc_treino,acc_teste=accuracy)) #retorna o primeiro mtry com a acurácia máxima
}

set.seed(123)
nbootstrap <- 10 #numero de reamostras # I used up to 10000 nbootstraps!
bootstrap_size <- 10 #genes por conjunto # I used from 10 to 50 bootstrap sizes!
indices<-lapply(1:nbootstrap,function(b)sample(1:nrow(data),bootstrap_size,replace=TRUE))
names(indices)<- paste0("Bootstrap",1:nbootstrap)
indices2<-lapply(indices, function(boot){
  as.integer(c(boot,boot+nrow(data)))#,boot+(nrow(data)*2)
})
value_cols <- 2:5
cols <- as.data.frame(t(combn(value_cols,3)))

accuracy_df <- data.frame()
models_boot <- list()
for(i in seq_along(indices2)) {
  cols_df <- data.frame()
  models_boot[[paste0("boot",i)]] <- list()
  for(c in 1:nrow(cols)) {
    boot_result <- calculate_accuracy(data=B4, indices=indices2[[i]], cols=as.numeric(cols[c,]))
    models_boot[[paste0("boot",i)]][[c]] <- boot_result$model
    cols_df <- rbind(cols_df,
                     cbind(boot=i,cols=c,acc_treino=boot_result$acc_treino,
                           acc_teste=boot_result$acc_teste))
  }
  cols_df$mean_acc_test <- median(cols_df$acc_teste)
  accuracy_df <- rbind(accuracy_df,cols_df)
}

#--------------Results from Classification are in the directory, including:
  #accuracy_df_10_10000 to accuracy_df_50_10000 , with the overview of all the results!
  #cols , which are the columns(samples) used for testing!
  #indices_10_10000 to indices_50_10000, which are the genes selected for each bootstrap!
 

#------Plot to see the results from classification: -------
acc10_B <- data_frame(Resample=accuracy_df_10_10000$boot,
                      Accuracy=accuracy_df_10_10000$mean_acc_test)
acc20_B <- data_frame(Resample=accuracy_df_20_10000$boot,
                      Accuracy=accuracy_df_20_10000$mean_acc_test)
acc30_B <- data_frame(Resample=accuracy_df_30_10000$boot,
                      Accuracy=accuracy_df_30_10000$mean_acc_test)
acc40_B <- data_frame(Resample=accuracy_df_40_10000$boot,
                      Accuracy=accuracy_df_40_10000$mean_acc_test)
acc50_B <- data_frame(Resample=accuracy_df_50_10000$boot,
                      Accuracy=accuracy_df_50_10000$mean_acc_test)
#acc100_B <- data_frame(Resample=accuracy_df_100_1000$boot,
#                     Accuracy=accuracy_df_100_1000$mean_acc_test)

detach(package:plyr)
acc10_B <- acc10_B %>% group_by(Resample) %>% 
  summarise(Accuracy = mean(Accuracy)) %>%
  filter(Accuracy >= 0.5) %>%
  mutate(Size=10)
acc20_B <- acc20_B %>% group_by(Resample) %>% 
  summarise(Accuracy = mean(Accuracy)) %>%
  filter(Accuracy >= 0.5) %>%
  mutate(Size=20)
acc30_B <- acc30_B %>% group_by(Resample) %>% 
  summarise(Accuracy = mean(Accuracy)) %>%
  filter(Accuracy >= 0.5) %>%
  mutate(Size=30)
acc40_B <- acc40_B %>% group_by(Resample) %>% 
  summarise(Accuracy = mean(Accuracy)) %>%
  filter(Accuracy >= 0.5) %>%
  mutate(Size=40)
acc50_B <- acc50_B %>% group_by(Resample) %>% 
  summarise(Accuracy = mean(Accuracy)) %>%
  filter(Accuracy >= 0.5) %>%
  mutate(Size=50)

acc_B <- rbind(acc10_B,acc20_B,acc30_B,acc40_B,acc50_B)

ggplot(data=acc_B, aes(x=Resample, y=Accuracy)) +
  facet_wrap(~Size, nrow = 1) +
  stat_identity(geom="line", position="dodge", colour="#0073e6") +
  coord_cartesian(ylim=c(0.5, 0.8)) +
  geom_hline(yintercept=0.7, linetype="dashed", color = "#f57600", size=1) +
  ylab("Median Test Accuracy") +
  theme_bw()
#WITH MEDIANS AND ACCURACIES FOR EACH SAMPLE TESTED  
acc_10_B_col <- accuracy_df_10_10000 %>% mutate(Size=10)
acc_20_B_col <- accuracy_df_20_10000 %>% mutate(Size=20)
acc_30_B_col <- accuracy_df_30_10000 %>% mutate(Size=30)
acc_40_B_col <- accuracy_df_40_10000 %>% mutate(Size=40)
acc_50_B_col <- accuracy_df_50_10000 %>% mutate(Size=50)
acc_B_col <- rbind(acc_10_B_col,acc_20_B_col,acc_30_B_col,acc_40_B_col,acc_50_B_col)
acc_B_col$cols <- cut(acc_B_col$cols, breaks = 4, labels = c("Sample 4", "Sample 3", "Sample 2", "Sample 1"))
acc_10_B_col$cols <- cut(acc_10_B_col$cols, breaks = 4, labels = c("Sample 4", "Sample 3", "Sample 2", "Sample 1"))

# Calculate the minimum and maximum acc_teste for each boot
acc_B_range <- acc_B_col %>% 
  group_by(boot, Size) %>% 
  summarize(acc_min = min(acc_teste), acc_max = max(acc_teste))

ggplot(acc_B_col) + 
  facet_wrap(~Size, nrow = 1) +
  geom_segment(data = acc_B_range, aes(x=boot, xend=boot, y=acc_min, yend=acc_max), color="lightgrey") +
  geom_point(aes(y = acc_teste, x=boot, color = cols), position="dodge") +
  geom_hline(yintercept=0.7, linetype="dashed", color = "black", size=0.5) +
  stat_identity(aes(y = mean_acc_test, x=boot), geom="line", position="dodge", color = "black", alpha=0.5) + #size=0.5, shape="triangle"
  coord_cartesian(ylim=c(0.1, 1)) +
  #coord_cartesian(xlim=c(1, 10)) +
  theme_bw() +
  theme(legend.position="bottom") +
  scale_color_manual(name = "Test", values = c("#0073e6", "#5ba300", "#f57600", "purple")) +
  xlab("") +
  ylab("Accuracy")

#FOCUSING ONLY IN THE RESAMPLES with bigger accuracies ---

acc_B_col_chosen <- acc_B_col %>%
  group_by(boot, Size) %>%
  filter(any(acc_teste >= 0.7))
# Calculate the minimum and maximum acc_teste for each boot
acc_B_range_chosen <- acc_B_col_chosen %>% 
  group_by(boot, Size) %>% 
  summarize(acc_min = min(acc_teste), acc_max = max(acc_teste))

ggplot(acc_B_col_chosen) + 
  facet_wrap(~Size, nrow = 1) +
  geom_segment(data = acc_B_range_chosen, aes(x=boot, xend=boot, y=acc_min, yend=acc_max), color="lightgrey") +
  geom_point(aes(y = acc_teste, x=boot, color = cols), position="dodge") +
  geom_hline(yintercept=0.7, linetype="dashed", color = "black", size=0.5) +
  stat_identity(aes(y = mean_acc_test, x=boot), geom="line", position="dodge", color = "black", alpha=0.5) + #size=0.5, shape="triangle"
  coord_cartesian(ylim=c(0.1, 1)) +
  #coord_cartesian(xlim=c(1, 10)) +
  scale_y_continuous(limits = c(0.1, 1), breaks = seq(0.1, 1, 0.1)) +
  theme_bw() +
  theme(legend.position="bottom") +
  scale_color_manual(name = "Test", values = c("#0073e6", "#5ba300", "#f57600", "purple")) +
  xlab("Resample") +
  ylab("Accuracy")

#FOCUSING ONLY IN THE 6 RESAMPLES CHOSEN ---
acc_B_col_chosen_6 <- acc_B_col %>%
  group_by(boot, Size) %>%
  filter(any(mean_acc_test >= 0.7))
# Calculate the minimum and maximum acc_teste for each boot
acc_B_range_chosen_6 <- acc_B_col_chosen_6 %>% 
  group_by(boot, Size) %>% 
  summarize(acc_min = min(acc_teste), acc_max = max(acc_teste))

ggplot(acc_B_col_chosen_6) + 
  facet_wrap(~Size, nrow = 1) +
  geom_segment(data = acc_B_range_chosen_6, aes(x=boot, xend=boot, y=acc_min, yend=acc_max), color="lightgrey") +
  geom_point(aes(y = acc_teste, x=boot, color = cols), position="dodge") +
  #geom_hline(yintercept=0.7, linetype="dashed", color = "black", size=0.5) +
  stat_identity(aes(y = mean_acc_test, x=boot), geom="line", position="dodge", color = "black", alpha=0.5) + #size=0.5, shape="triangle"
  coord_cartesian(ylim=c(0.1, 1)) +
  #coord_cartesian(xlim=c(1, 10)) +
  scale_y_continuous(limits = c(0.1, 1), breaks = seq(0.1, 1, 0.1)) +
  theme_bw() +
  theme(legend.position="bottom") +
  scale_color_manual(name = "Test", values = c("#0073e6", "#5ba300", "#f57600", "purple")) +
  xlab("Resample") +
  ylab("Accuracy")

#Get the BEST GENES after classification: ----------------
# BONE 4 2C: Getting the best GENES (UNIQUE or UNION) of the best boots from the classification:
#Best boots: 166,799,2742,5019,5449,7064 <- They all showed 0.7 median test accuracy!
#50 genes unique!
melhor_genes_unique_B <- cbind(indices_10_10000[[166]],indices_10_10000[[799]],indices_10_10000[[2742]],indices_10_10000[[5019]],indices_10_10000[[5449]],indices_10_10000[[7064]])
melhor_genes_unique_B <- as.data.frame(melhor_genes_unique_B)
melhor_genes_unique_B <- unique(c(melhor_genes_unique_B[,1],melhor_genes_unique_B[,2], melhor_genes_unique_B[,3], melhor_genes_unique_B[,4], melhor_genes_unique_B[,5], melhor_genes_unique_B[,6]))


melhor_genes_B <- B4[melhor_genes_unique_B,]
melhor_genes_B <- melhor_genes_B[[1]]
melhor_genes_B <- top_genes_B4[,melhor_genes_B]
pca_melhor_B4_2C <- prcomp(melhor_genes_B, center = TRUE, scale = FALSE) # center = TRUE ,scale = FALSE
summary(pca_melhor_B4_2C)

#OPTION: Choosing the bootstraps with the smaller SD: 
#Getting maybe even better bootstraps among those 6 already better... then I will get 
#the bootstraps among those 6 that have the smaller SD:
#Acc_tests obtained with each boot:
# sd(accuracy_df_10_10000$acc_teste[661:664]) #Boot 166
# [1] 0.1108678
# > sd(accuracy_df_10_10000$acc_teste[3193:3196])
# [1] 0.1887459
# > sd(accuracy_df_10_10000$acc_teste[10965:10968])
# [1] 0.2015564
# > sd(accuracy_df_10_10000$acc_teste[20073:20076])
# [1] 0.1108678
# > sd(accuracy_df_10_10000$acc_teste[21793:21796])
# [1] 0.1108678
# > sd(accuracy_df_10_10000$acc_teste[28253:28256])
# [1] 0.1414214

#--Smaller SD: 166,5019,5449---
#27 genes!
melhor_genes_unique_B <- cbind(indices_10_10000[[166]],indices_10_10000[[5019]],indices_10_10000[[5449]])
melhor_genes_unique_B <- as.data.frame(melhor_genes_unique_B)
melhor_genes_unique_B <- unique(c(melhor_genes_unique_B[,1],melhor_genes_unique_B[,2], melhor_genes_unique_B[,3]))
melhor_genes_unique_B <- intersect(melhor_genes_unique_B[,2],melhor_genes_unique_B[,3])

melhor_genes_B <- B4[melhor_genes_unique_B,]
melhor_genes_B <- melhor_genes_B[[1]]
melhor_genes_B <- top_genes_B4[,melhor_genes_B]
#melhor_genes_B <- colnames(melhor_genes_B)
pca_melhor_B4 <- prcomp(melhor_genes_B, center = TRUE, scale = FALSE) # center = TRUE ,scale = FALSE
summary(pca_melhor_B4)

#OPTION: Getting the worst genes:
#WORST GENES
#8122 and 9532
piores_genes_unique_B <- cbind(indices_10_10000[[8122]],indices_10_10000[[9532]])
piores_genes_unique_B <- as.data.frame(piores_genes_unique_B)
piores_genes_unique_B <- unique(c(piores_genes_unique_B[,1],piores_genes_unique_B[,2]))


pior_genes_B <- B4[piores_genes_unique_B,]
pior_genes_B <- pior_genes_B[[1]]
pior_genes_B <- top_genes_B4[,pior_genes_B]
pca_pior_B4 <- prcomp(pior_genes_B, center = TRUE, scale = FALSE) # center = TRUE ,scale = FALSE
summary(pca_pior_B4)

#---- Getting gene names instead of Gene IDs -------------------
#Using the gene_IDS or the gene_names using the GTF data for the organism to get the gene_names
org.Mm.eg.db #To get the details of this annotation

#GTF_data <- unique(GTF_data) #GTF_data from Deseq2 BoneSkin3.R
#melhor_genes_10_10000 <- colnames(melhor_genes_10_10000)
melhor_genes_B <- colnames(melhor_genes_B)
#melhor_genes_50_10000 <- colnames(melhor_genes_50_10000)

melhor_genes_B4 <- sub("\\.\\d+$", "", melhor_genes_B) #melhor_genes_B from Supervised.R
#melhor_genes_50_B4 <- sub("\\.\\d+$", "", melhor_genes_50_10000) #Melhor_genes_10_10000 from Supervised.R
#pior_genes_10_B4 <- sub("\\.\\d+$", "", pior_genes_10_10000) #Melhor_genes_10_10000 from Supervised.R

##------ Preparing the GTF data to get the names for the genes_ids
gtf_data <- gencode.vM35.primary_assembly.annotation #THis is the gtf data uploaded!
# Split the ninth column based on semicolon (;)
gtf_data$split_data <- strsplit(gtf_data[, 9], ";")

# Separate gene_id and gene_name from the split data
gtf_data$gene_id <- sapply(gtf_data$split_data, function(x) x[grepl("gene_id", x)])
gtf_data$gene_name <- sapply(sapply(gtf_data$split_data, function(x) x[grepl("gene_name", x)]), function(x) str_replace_all(x, "gene_name ", ""))

# Clean gene_name (remove extra spaces)
gtf_data$gene_name <- gsub(" +$", "", gtf_data$gene_name)

# Remove the temporary split column
gtf_data <- gtf_data[ , !(names(gtf_data) %in% "split_data")]
GTF_data <- gtf_data[,10:11]
# Remove "gene_id " from the ninth column
GTF_data[, 1] <- gsub("^gene_id ", "", GTF_data[, 1])

#To use later, getting the names of the classification genes:
GTF_data <- unique(GTF_data)
melhor_genes_B_gtf <- GTF_data %>% filter(GTF_data$gene_id %in% melhor_genes_B)

#---- Comparison between tissues. If there is more than 1 tissue, a GENE intersection can be made:-------
#Getting the genes that are used for classification for both tissues that are the same or just to that tissue:
Same_Genes_Tissues <- intersect(melhor_genes_S_gtf$gene_name, melhor_genes_B_gtf$gene_name)
intersect(melhor_genes_S_gtf$gene_id, melhor_genes_B_gtf$gene_id)

#Getting the unique genes (the difference between the first object and the second object)
Unique_Genes_S4 <- setdiff(melhor_genes_S_gtf$gene_name, melhor_genes_B_gtf$gene_name) #For unique in S
Unique_Genes_B4 <- setdiff(melhor_genes_B_gtf$gene_name, melhor_genes_S_gtf$gene_name)  # For unique in B

#After classification... LOG2FC of the selected GENES -------
#Using DESeq2 package....
log2FC_B_SxF <- as.data.frame(res_B4_SxF)
log2FC_B_SxF_filtered <- log2FC_B_SxF %>% filter(log2FoldChange >= 1.5 | log2FoldChange <= -1.5)
log2FC_B_SxF <- log2FC_B_SxF[, c("log2FoldChange", "lfcSE")]

#Using the genes selected after classification
melhor_genes_B <- melhor_genes_B_gtf %>% arrange((gene_id))

log2fc_B <- log2FC_B_SxF %>%
  filter(row.names(log2FC_B_SxF) %in% melhor_genes_B$gene_id) %>%
  mutate(mycolor= ifelse(log2FoldChange>0, "#positive", "#negative"))
rownames_B_log2fc <- row.names(log2fc_B)
#rownames_B_log2fc <- sub("\\.\\d+$", "", rownames_B_log2fc) #melhor_genes_B from Supervised.R
log2fc_B <- log2fc_B %>%
  mutate(gene_id = rownames_B_log2fc) %>%
  arrange((gene_id)) %>%
  mutate(gene_name = melhor_genes_B$gene_name)


ggplot(log2fc_B) +
  geom_hline(yintercept=0, linetype="solid", color = "black", size=0.2) +
  geom_segment(aes(x=gene_name, xend=gene_name, y=0, yend=log2FoldChange), color="black") +
  geom_point(aes(y = log2FoldChange, x=gene_name, color = mycolor), position="dodge", size = 3) +
  #  stat_identity(aes(x =gene_name, y = log2FoldChange, color = mycolor))+
  coord_flip() +
  scale_y_continuous(limits = c(-3, 3), breaks = seq(-3, 3, 0.5)) +
  scale_color_manual(values = c("#f57600", "#0073e6")) +
  xlab("") +
  ylab("Log2FC")+
  theme_bw() +
  theme(legend.position = "none")

#Correlation analysis-----------------
corr_10_B4 <- top_genes_B4[,melhor_genes_B$gene_id]
corr_10_B4 <- as.data.frame(corr_10_B4)

corr_10_B4_F <- corr_10_B4[1:4,]
corr_10_B4_S <- corr_10_B4[5:8,]

#Changing the gene_ids to gene_names
corr_10_B4_name <- corr_10_B4 %>%
  rename_all(~ melhor_genes_B_gtf$gene_name[match(., melhor_genes_B_gtf$gene_id)])

corr_10_B4_name_F <- corr_10_B4_name[1:4,]
corr_10_B4_name_S <- corr_10_B4_name[5:8,]

#-----Correlation between same genes from different conditions (F and S):
# Calculate correlations for each gene
correlations_10_B4 <- sapply(1:ncol(corr_10_B4), function(i) cor(corr_10_B4_F[, i], corr_10_B4_S[, i], method = "spearman"))
gene_names_corr_10_B4 <- colnames(corr_10_B4_name)
gene_pairs_corr_10_B4 <- cbind(gene_names_corr_10_B4, correlations_10_B4)
colnames(gene_pairs_corr_10_B4) <- c("Gene", "Correlation")
gene_pairs_corr_10_B4 <- as.data.frame(gene_pairs_corr_10_B4)
gene_pairs_corr_10_B4$Correlation <- as.numeric(gene_pairs_corr_10_B4$Correlation)

# Barplot
gene_pairs_corr_10_B4 <- gene_pairs_corr_10_B4 %>%
  mutate(mycolor= ifelse(Correlation>0, "#0073e6", "#f57600")) %>%
  arrange(Correlation)  %>%   # First sort by val. This sort the dataframe but NOT the factor levels
  mutate(Gene=factor(Gene, levels=Gene))# This trick update the factor levels

ggplot(gene_pairs_corr_10_B4, aes(x=Gene, y=Correlation)) + 
  geom_bar(stat = "identity", position = "dodge", colour = gene_pairs_corr_10_B4$mycolor, fill = gene_pairs_corr_10_B4$mycolor, width = 0.7) +
  coord_flip() +
  #Above 0.7 starts a very strong relationship!!!!
  geom_hline(yintercept=0.7, linetype="dashed", color = "#0073e6", size=1) +
  geom_hline(yintercept=-0.7, linetype="dashed", color = "#f57600", size=1) +
  theme_bw()

#------------ Correlation between GENES for each condition
correlation_matrix_10_B4 <- cor(corr_10_B4) ##NOT USED: Correlation between genes but using both conditions together!

#now focusing in each condition for comparison
correlation_matrix_10_B4_F <- cor(corr_10_B4_F)
colnames(correlation_matrix_10_B4_F) <- colnames(corr_10_B4_name)
rownames(correlation_matrix_10_B4_F) <- colnames(corr_10_B4_name)
# Create a correlation heatmap with order parameter using known AOE method
corrplot(correlation_matrix_10_B4_F, method = "color", col=colorRampPalette(c("#f57600","white","#0073e6"))(100), type = "lower", tl.cex = 0.5, tl.col="black", number.cex = 0.5, order = "AOE") #To add values inside: addCoef.col = "black", #For significance: p.mat = cor_test_S4_F$p, sig.level = 0.05,

correlation_matrix_10_B4_S <- cor(corr_10_B4_S)
colnames(correlation_matrix_10_B4_S) <- colnames(corr_10_B4_name)
rownames(correlation_matrix_10_B4_S) <- colnames(corr_10_B4_name)
# Create a correlation heatmap with order parameter using known AOE method
corrplot(correlation_matrix_10_B4_S, method = "color", col=colorRampPalette(c("#f57600","white","#0073e6"))(100), type = "lower", tl.cex = 0.5, tl.col="black", number.cex = 0.5, order = "AOE") #To add values inside: addCoef.col = "black", #For significance: p.mat = cor_test_S4_F$p, sig.level = 0.05,

# NOT USED in this research: Top most correlated values -- But it uses pvalue!
install.packages("lares")
library(lares)
corr_cross(corr_10_B4, rm.na = T, max_pvalue = 0.05, top = 10, grid = T)
corr_cross(corr_10_B4_F, rm.na = T, max_pvalue = 0.05, top = 10, grid = T)
corr_cross(corr_10_B4_S, rm.na = T, max_pvalue = 0.05, top = 10, grid = T)


#Getting the intersect and UNIQUE correlations that are used for conditions in each tissue:
cor_diff_B4 <- abs(correlation_matrix_10_B4_F - correlation_matrix_10_B4_S)

#Using a cutoff of the absolute difference between correlations
#cor_diff_unique_B4 <- cor_diff_B4 >= 0.3 #Choose the cutoff!!!!
#cor_diff_B4[!cor_diff_unique_B4] <- NA

#plotting
pheatmap(cor_diff_B4, breaks = c(0,0.5,1,1.5,2), col = colorRampPalette(c("#fffd8d", "yellow", "red"))(4), na.rm = TRUE)
#dev.off()

#Creating a table with the correlations with differences bigger than 0.3
# Get the indices of non-NA values in filtered_cor_diff
indices_cor_B4 <- which(!is.na(cor_diff_B4), arr.ind = TRUE)

# Create a data frame with gene names and correlation differences
result_table_cor_B4 <- data.frame(
  Gene1 = rownames(cor_diff_B4)[indices_cor_B4[, 1]],
  Gene2 = colnames(cor_diff_B4)[indices_cor_B4[, 2]],
  Difference = cor_diff_B4[indices_cor_B4]
)

#plot alternative
ggplot(result_table_cor_B4, aes(x=Gene1, y=Gene2, fill=Difference)) +
  geom_tile() +
  scale_fill_gradient2(low = "#fffd8d", mid = "yellow", high = "red") +
  #scale_fill_gradientn(colors = hcl.colors(20, "YlOrRd")) +
  theme_bw() +
  theme(legend.position = "bottom")


#------------Association Rule Mining: -----------------
# FOR FLASH:
# Assuming 'counts_matrix' is your 10x4 matrix
transactions_B4_F <- as(corr_10_B4_name_F, "transactions") #Use for gene names: corr_10_B4_name_F
dim(transactions_B4_F) #We have A transactions (baskets) x B items in each basket
itemLabels(transactions_B4_F) #The genes in each basket
#max(corr_10_B4_F$ENSMUSG00000061315.15)
#min(corr_10_B4_F$ENSMUSG00000061315.15)

summary(transactions_B4_F)
#Showing the ocorrence of the genes in each sample:
#We can see the level of each gene, each sample presents:
image(transactions_B4_F)
#Display the relative item frequency:
itemFrequencyPlot(transactions_B4_F, topN=10,  cex.names=0.5, mai = c(2, 1, 0.5, 0))
#sort(itemFrequency(transactions_B4_F, type="relative")) 


# Mine association rules with minimum support of 0.5 and minimum confidence of 0.8
#With small data (increase the support(frequency) that the data appear to ensure higher confidence)):
#Both apriori and FP-Growtn method gave the same results! FP-Growth could deal better with variation, 
#however it makes the difference in larger data sizes!
rules_B4_F <- apriori(transactions_B4_F, parameter = list(supp=0.5, conf=0.8, minlen=2, maxlen=10, target= "rules"))#, algorithm = "fpgrowth"
#fpgrowth_B4_F <- arules::fim4r(transactions_B4_F, method = "fpgrowth", target = "rules", supp = 0.5, conf = 0.8, zmin = 2, zmax = 10)
summary(rules_B4_F)
#summary(fpgrowth_B4_F)

# Rules with 1 gene were ignored (minlen=2), since:
#   * an empty LHS mean that no matter what other items are involved the item in the RHS will appear with the probability given by the rule’s confidence (which equals the support) *

#Too many rules, need to filter: 
#By lift, however, all my rules had the same lift = 2 (good)! 
##select(head(sort(fpgrowth_B4_F, by ="lift"),10))

#Removing redundant rules: 
# a rule is considered redundant if it is equally or less predictive than a more general rule,
#which has the same items on the right hand side, 
#but one or more items less in the left hand side. 
rules_B4_F <- rules_B4_F[!is.redundant(rules_B4_F)]   #from almost 700 thousand to 500 rules
#fpgrowth_B4_F <- fpgrowth_B4_F[!is.redundant(fpgrowth_B4_F)]   #from almost 800 thousand to 500 rules

#Still a lot of rules: 
#Removing statistically insignificant rules:
# NOT USED, but when tested yield the same results. 
rules_B4_F <- rules_B4_F[!is.significant(rules_B4_F, 
                                         transactions_B4_F, 
                                         method = "fisher", 
                                         adjust = 'bonferroni')]
# fpgrowth_B4_F <- fpgrowth_B4_F[!is.significant(fpgrowth_B4_F, 
#                                          transactions_B4_F, 
#                                          method = "fisher", 
#                                          adjust = 'bonferroni')]
summary(rules_B4_F)
#summary(fpgrowth_B4_F)

#If needed: Focus on rules for an item in particular:
# ## Extract rules that have "Churn" as consequent 
# churn_rules <- subset(sig_rules, subset=rhs %pin% 'Churn')
# 
# summary(churn_rules)

# Inspect the rules
inspect(rules_B4_F)
#inspect(fpgrowth_B4_F)

#Turning in to dataframes:
rules_B4_F_df <- DATAFRAME(rules_B4_F, setStart='', setEnd='', separate = TRUE)
#fpgrowth_B4_F_df <- DATAFRAME(fpgrowth_B4_F, setStart='', setEnd='', separate = TRUE)

#Ways to visualize large set of rules: Scatter plots, Graph and grouped matrix
# https://cran.r-project.org/web/packages/arulesViz/vignettes/arulesViz.pdf 

#Function plot of the package arulezviz
#Scatterplot, however for no variation in parameters it doesnt help us
plot(rules_B4_F)
#Interactive plot of the rules
plot(rules_B4_F, engine = "plotly")

#Showing all the rules for that condition in hubs: INTERACTIVE, nice!
#subrules_B4_F <- head(rules_B4_F, n = 596, by = "confidence")
plot(rules_B4_F, method = "graph",  engine = "htmlwidget", max = 600)

#This graph with igraph plot is not good for too many rules!
#plot(rules_B4_F, method = "graph",  engine = "igraph", max = 600) #, layout = igraph::in_circle()

#This paracoord is not good for too many rules!
#plot(rules_B4_F, method="paracoord", max = 600)

#Nice options for plotting per condition: Showing ALL THE RULES using HUBS!!! How many rules have that LHS/gene in that level:
plot(rules_B4_F, method = "grouped matrix", measure = "support", shading = "lift", control = list(k = 50), rhs_max	 =  50, main = NULL, col = "#0073e6") 
plot(rules_B4_F, method = "matrix", reorder = "measure")

#Showing ALL rules: Visualizing with Gephi (installed on my PC (for all users)):
#saveAsGraph(head(rules_B4_F, n = 596, by = "lift"), file = "rules.graphml")
#library("igraph")
#saveAsGraph(rules_B4_F, file = "rules_B4_F.graphml")

# STANDARD:

# Assuming 'counts_matrix' is your 10x4 matrix
transactions_B4_S <- as(corr_10_B4_name_S, "transactions")
dim(transactions_B4_S) #We have A transactions (baskets) x B items in each basket
itemLabels(transactions_B4_S) #The genes in each basket
#max(corr_10_B4_F$ENSMUSG00000061315.15)
#min(corr_10_B4_F$ENSMUSG00000061315.15)

summary(transactions_B4_S)
#Showing the ocorrence of the genes in each sample:
#We can see the level of each gene, each sample presents:
image(transactions_B4_S)
#Display the relative item frequency:
itemFrequencyPlot(transactions_B4_S, topN=10,  cex.names=0.45, mai = c(1, 1, 1, 1))


# Mine association rules with minimum support of 0.5 and minimum confidence of 0.8
rules_B4_S <- apriori(transactions_B4_S, parameter = list(supp=0.5, conf=0.8,  minlen=2, maxlen=10, target= "rules"))
#fpgrowth_B4_S <- fim4r(transactions_B4_S, method = "fpgrowth", target = "rules", supp = 10, conf = 10)

summary(rules_B4_S)


#Removing redundant rules: 
#lift was not neccessary, as were all the same. Lift =2 (good)!

#Using is.redundant function for removal of redundant rules:
rules_B4_S <- rules_B4_S[!is.redundant(rules_B4_S)]   #from almost 700 thousand to 500 rules
#fpgrowth_B4_F <- fpgrowth_B4_F[!is.redundant(fpgrowth_B4_F)]   #from almost 800 thousand to 500 rules

#Still a lot of rules: 
#NOT USED, but yielded no difference in the results when used in this case: Removing statistically insignificant rules:
rules_B4_S <- rules_B4_S[!is.significant(rules_B4_S, 
                                         transactions_B4_S, 
                                         method = "fisher", 
                                         adjust = 'bonferroni')]

summary(rules_B4_S)

# Inspect the rules
inspect(rules_B4_S)
#inspect(fpgrowth_B4_F)

#Turning in to dataframes:
rules_B4_S_df <- DATAFRAME(rules_B4_S, setStart='', setEnd='', separate = TRUE)
#fpgrowth_B4_F_df <- DATAFRAME(fpgrowth_B4_F, setStart='', setEnd='', separate = TRUE)

#PLotting but using lift. In our case its not neccessary:
plot(rules_B4_S)
#Interactive plot of the rules. INTERACTIVE!
plot(rules_B4_S, engine = "plotly")

#subrules_B4_S <- head(rules_B4_S, n = 476, by = "confidence")

#plot(subrules_B4_S, method = "graph",  engine = "htmlwidget")

#plot(subrules_B4_S, method="paracoord")


#Showing all the rules per condition in hubs, INTERACTIVE, nice!!!
#subrules_B4_F <- head(rules_B4_F, n = 596, by = "confidence")
plot(rules_B4_S, method = "graph",  engine = "htmlwidget", max = 600)

#This graph with igraph plot is not good for too many rules!
#plot(rules_B4_F, method = "graph",  engine = "igraph", max = 600) #, layout = igraph::in_circle()

#This paracoord is not good for too many rules!
#plot(rules_B4_F, method="paracoord", max = 600)

#Showing ALL THE RULES using HUBS.NICE!!! How many rules have that LHS/gene in that level:
plot(rules_B4_S, method = "grouped matrix", measure = "support", shading = "lift", control = list(k = 50), rhs_max	 =  50, main = NULL, col = "#0073e6") 
plot(rules_B4_S, method = "matrix", reorder = "measure")

#Showing ALL rules: Visualizing with Gephi (installed on my PC (for all users)):
#saveAsGraph(head(rules_B4_F, n = 596, by = "lift"), file = "rules.graphml")
#library("igraph")
#saveAsGraph(rules_B4_S, file = "rules_B4_S.graphml")

#Unique and Shared RULES: Comparison of rules between Conditions (Shared and Unique for each condition) -----
gene_rules_B4_F <- rules_B4_F_df
gene_rules_B4_F$LHS <- gsub("\\=.*", "", gene_rules_B4_F$LHS)
gene_rules_B4_F$RHS <- gsub("\\=.*", "", gene_rules_B4_F$RHS)

gene_rules_B4_S <- rules_B4_S_df
gene_rules_B4_S$LHS <- gsub("\\=.*", "", gene_rules_B4_S$LHS)
gene_rules_B4_S$RHS <- gsub("\\=.*", "", gene_rules_B4_S$RHS)

intersect_gene_rules_B4 <- intersect(gene_rules_B4_F, gene_rules_B4_S)
unique_gene_rules_B4_F <- setdiff(gene_rules_B4_F, gene_rules_B4_S) #For FLASH only
unique_gene_rules_B4_S <- setdiff(gene_rules_B4_S, gene_rules_B4_F) #For STANDARD only

intersect_gene_rules_B4 <- intersect_gene_rules_B4[,1:2]
unique_gene_rules_B4_F <- unique_gene_rules_B4_F[,1:2]
unique_gene_rules_B4_S <- unique_gene_rules_B4_S[,1:2]

# -------- Create a graph from the data frame (INTERSECT)
library("igraph")
graph_gene_rules_B4 <- graph_from_data_frame(intersect_gene_rules_B4, directed = TRUE)

# Get the cluster memberships
com_B4_intersect <- cluster_walktrap(graph_gene_rules_B4)
membership <- membership(com_B4_intersect)

# Add the cluster membership as a new column to the data frame
intersect_gene_rules_B4$grp <- membership[V(graph_gene_rules_B4)$name]

# Get unique group IDs
unique_groups_intersect_B4 <- unique(intersect_gene_rules_B4$grp)

# Create a color palette
mycolor <- brewer.pal(n = length(unique_groups_intersect_B4), name = "Paired")

# Assign colors to nodes based on their group membership
V(graph_gene_rules_B4)$color <- mycolor[membership]

#Arc diagram
library("ggraph")
ggraph(graph_gene_rules_B4, layout="linear") +
  geom_edge_arc(edge_colour="black", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
  geom_node_point(aes(size=5, color=as.factor(color), fill=color), alpha=0.5) +
  scale_size_continuous(range=c(0.5,8)) +
  scale_color_manual(values=mycolor) +
  geom_node_text(aes(label=name), angle=65, hjust=1, nudge_y = -0.3, size=2.3) +
  theme_void() +
  theme(legend.position="none",
        plot.margin=unit(c(0,0,0.4,0), "null"),
        panel.spacing=unit(c(0,0,3.4,0), "null")) +
  expand_limits(x = c(-1.2, 1.2), y = c(-5.6, 1.2))


#------------ Create a graph from the data frame (Unique for FLASH)
graph_unique_gene_rules_B4_F <- graph_from_data_frame(unique_gene_rules_B4_F, directed = TRUE)

# Get the cluster memberships
com_B4_unique_F <- cluster_walktrap(graph_unique_gene_rules_B4_F)
membership_B4_F <- membership(com_B4_unique_F)

membership_B4_F <- membership_B4_F[order(membership_B4_F)]

#Reorder dataset and make the graph
unique_gene_rules_B4_F <- unique_gene_rules_B4_F %>%
  arrange(match(LHS, names(membership_B4_F))) %>%
  mutate(LHS=factor(LHS))

# Create a graph from the data frame (Unique for FLASH)
graph_unique_gene_rules_B4_F <- graph_from_data_frame(unique_gene_rules_B4_F, directed = TRUE)

# Get the cluster memberships
com_B4_unique_F <- cluster_walktrap(graph_unique_gene_rules_B4_F)
membership_B4_F <- membership(com_B4_unique_F)

# Add the cluster membership as a new column to the data frame
#unique_gene_rules_B4_F$grp <- membership_B4_F[V(graph_unique_gene_rules_B4_F)$name]

# Get unique group IDs
unique_groups_unique_B4_F <- unique(membership_B4_F)

# Create a color palette
mycolor <- brewer.pal(n = length(unique_groups_unique_B4_F), name = "Paired")

# Assign colors to nodes based on their group membership
V(graph_unique_gene_rules_B4_F)$color <- mycolor[membership_B4_F]



#Arc diagram
ggraph(graph_unique_gene_rules_B4_F, layout="linear") +
  geom_edge_arc(edge_colour="black", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
  geom_node_point(aes(size=5, color=as.factor(color), fill=color), alpha=0.5) +
  scale_size_continuous(range=c(0.5,8)) +
  scale_color_manual(values=mycolor) +
  geom_node_text(aes(label=name), angle=65, hjust=1, nudge_y = -0.3, size=2.3) +
  theme_void() +
  theme(legend.position="none",
        plot.margin=unit(c(0,0,0.4,0), "null"),
        panel.spacing=unit(c(0,0,3.4,0), "null")) +
  expand_limits(x = c(-1.2, 1.2), y = c(-5.6, 1.2))

# ----------- Create a graph from the data frame (Unique for STANDARD)
graph_unique_gene_rules_B4_S <- graph_from_data_frame(unique_gene_rules_B4_S, directed = TRUE)

# Get the cluster memberships
com_B4_unique_S <- cluster_walktrap(graph_unique_gene_rules_B4_S)
membership_B4_S <- membership(com_B4_unique_S)

membership_B4_S <- membership_B4_S[order(membership_B4_S)]

#Reorder dataset and make the graph
unique_gene_rules_B4_S <- unique_gene_rules_B4_S %>%
  arrange(match(LHS, names(membership_B4_S))) %>%
  mutate(LHS=factor(LHS))

# Create a graph from the data frame (Unique for FLASH)
graph_unique_gene_rules_B4_S <- graph_from_data_frame(unique_gene_rules_B4_S, directed = TRUE)

# Get the cluster memberships
com_B4_unique_S <- cluster_walktrap(graph_unique_gene_rules_B4_S)
membership_B4_S <- membership(com_B4_unique_S)

# Add the cluster membership as a new column to the data frame
#unique_gene_rules_B4_S$grp <- membership_B4_S[V(graph_unique_gene_rules_B4_S)$name]

# Get unique group IDs
unique_groups_unique_B4_S <- unique(membership_B4_S)

# Create a color palette
mycolor <- brewer.pal(n = length(unique_groups_unique_B4_S), name = "Paired")

# Assign colors to nodes based on their group membership
V(graph_unique_gene_rules_B4_S)$color <- mycolor[membership_B4_S]



#Arc diagram
ggraph(graph_unique_gene_rules_B4_S, layout="linear") +
  geom_edge_arc(edge_colour="black", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
  geom_node_point(aes(size=5, color=as.factor(color), fill=color), alpha=0.5) +
  scale_size_continuous(range=c(0.5,8)) +
  scale_color_manual(values=mycolor) +
  geom_node_text(aes(label=name), angle=65, hjust=1, nudge_y = -0.3, size=2.3) +
  theme_void() +
  theme(legend.position="none",
        plot.margin=unit(c(0,0,0.4,0), "null"),
        panel.spacing=unit(c(0,0,3.4,0), "null")) +
  expand_limits(x = c(-1.2, 1.2), y = c(-5.6, 1.2))




#------------Classification using different algorithms at once:----------
#Using the gene refined after classification: In this case there were 235 genes!
B4_F <- B4[1:235,1:6] #<----- (top 100), B4[1:117,1:6] (top 50),  B4[1:57,1:6] (top 25), B4[1:24,1:6] (top 10),  B4[1:2236,1:6](top 1000), ou B4_top_genes_imp[1:55,1:6] ou B4_all_genes[1:21590,1:6] ou B4[1:2236,1:6] (for top 1000 genes)
data <- B4_F

#----- Just 1 algorithm - Random Forest base (RF) ---


calculate_accuracy <- function(data, indices, cols, modelo) {
  train_data <- data[indices, c(1,cols,6)]
  test_data <- data[indices,!(1:length(data)) %in% c(cols,6)]
  test_data <- cbind(test_data,replicate(length(cols)-1,test_data[2]))
  test_data$Condition <- data[indices,"Condition"]
  colnames(test_data) <- c("Gene",names(data[cols]),"Condition")
  indices_train_data <- seq_along(indices) #todos os indices em train_data (bootstrapinho)
  set.seed(123)
  #  browser()
  model <- train(Condition~.,  
                 data=train_data,  
                 method=modelo,
                 metric = "Accuracy", 
                 preProcess = "center",
                 trControl=trainControl(method="boot", number=1,index = list(indices_train_data),indexOut = list(indices_train_data)),
                 tuneLength = 10)
  acc_treino <- model$results%>%
    filter(mtry==model$bestTune$mtry)%>%
    select(Accuracy)
  predictions <- predict(model, newdata = test_data[1:5])
  accuracy <- mean(predictions == test_data$Condition)
  return(list(model=model,acc_treino=acc_treino,acc_teste=accuracy)) #retorna o primeiro mtry com a acurácia máxima
}

set.seed(123)
nbootstrap <- 10 #numero de reamostras
bootstrap_size <- 10 #genes por conjunto
indices<-lapply(1:nbootstrap,function(b)sample(1:nrow(data),bootstrap_size,replace=TRUE))
names(indices)<- paste0("Bootstrap",1:nbootstrap)
indices2<-lapply(indices, function(boot){
  as.integer(c(boot,boot+nrow(data)))#,boot+(nrow(data)*2)
})
value_cols <- 2:5
cols <- as.data.frame(t(combn(value_cols,3)))

accuracy_df <- data.frame()
models_boot <- list()

#modelos <- names(getModelInfo())
modelos <- c("rf")
for(m in modelos) {
  for(i in seq_along(indices2)) {
    cols_df <- data.frame()
    models_boot[[paste0("boot",i)]] <- list()
    for(c in 1:nrow(cols)) {
      #boot_result <- tryCatch(calculate_accuracy(data=B4_top_genes, indices=indices2[[i]], cols=as.numeric(cols[c,]), modelo=m),error = function(e) NA)
      boot_result <- calculate_accuracy(data=B4_top_genes, indices=indices2[[i]], cols=as.numeric(cols[c,]), modelo=m)
      if(length(boot_result)<1) next
      models_boot[[paste0("boot",i)]][[c]] <- boot_result$model
      cols_df <- rbind(cols_df,
                       cbind(model=m,boot=i,cols=c,acc_treino=boot_result$acc_treino,
                             acc_teste=boot_result$acc_teste))
    }
    cols_df$mean_acc_test <- median(cols_df$acc_teste)
    accuracy_df <- rbind(accuracy_df,cols_df)
  }
  write.csv(accuracy_df,paste0("model_", m,".csv"))
}

---------------------- # Mais de um modelo principal para testar
# Create a list to store hyperparameters for each algorithm!!! <---- Add here!!!!
hyper_params <- list(
  rf = data.frame(mtry = c(1, 5, 10, sqrt(length(cols)))), 
  ranger = data.frame(mtry = c(1, 5, 10, sqrt(length(cols))), 
                      splitrule = "extratrees", 
                      min.node.size = c(2, 3)), 
  svmLinear2 = data.frame(cost = c(0.1, 1, 10)), 
  svmRadial = data.frame(sigma = c(0.1, 0.5, 1), C = c(0.1, 1, 10)), 
  nb = data.frame(fL = rep(0, 6), 
                  usekernel = rep(TRUE, 6), 
                  adjust = rep(1, 6)), 
  knn = data.frame(k = c(3, 5, 7)), 
  gbm = data.frame(interaction.depth = c(1, 2, 3), 
                   n.trees = c(50, 100, 150), 
                   shrinkage = c(0.01, 0.001, 0.005), 
                   n.minobsinnode = c(1, 2, 3)) 
)

calculate_accuracy <- function(data, indices, cols, modelo) {
  train_data <- data[indices, c(1,cols,6)]
  test_data <- data[indices,!(1:length(data)) %in% c(cols,6)]
  test_data <- cbind(test_data,replicate(length(cols)-1,test_data[2]))
  test_data$Condition <- data[indices,"Condition"]
  colnames(test_data) <- c("Gene",names(data[cols]),"Condition")
  indices_train_data <- seq_along(indices) #todos os indices em train_data (bootstrapinho)
  set.seed(123)
  #  browser()
  model <- train(Condition~.,  
                 data=train_data,  
                 method=modelo,
                 metric = "Accuracy", 
                 preProcess = "center",
                 trControl=trainControl(method="boot", number=1,index = list(indices_train_data),indexOut = list(indices_train_data)),
                 tuneGrid = hyper_params[[modelo]])
  acc_treino <- model$results%>%
    filter(row_number() == 1) %>% # Select the first row (usually the best performance)
    select(Accuracy)
  # acc_treino <- model$results%>%
  #   filter(mtry==model$bestTune$mtry)%>%
  #   select(Accuracy)
  predictions <- predict(model, newdata = test_data[1:5])
  accuracy <- mean(predictions == test_data$Condition)
  return(list(model=model,acc_treino=acc_treino,acc_teste=accuracy)) #retorna o primeiro mtry com a acurácia máxima
}

set.seed(123)
nbootstrap <- 10 #numero de reamostras # CHANGE HERE!
bootstrap_size <- 10 #genes por conjunto  # CHANGE HERE!
indices<-lapply(1:nbootstrap,function(b)sample(1:nrow(data),bootstrap_size,replace=TRUE))
names(indices)<- paste0("Bootstrap",1:nbootstrap)
indices2<-lapply(indices, function(boot){
  as.integer(c(boot,boot+nrow(data)))#,boot+(nrow(data)*2)
})
value_cols <- 2:5
cols <- as.data.frame(t(combn(value_cols,3)))

accuracy_df <- data.frame()
models_boot <- list()

#modelos <- names(getModelInfo())
#Choose algorithms here!!!!!!!!! :
modelos <- c("rf", "ranger", "svmLinear2", "svmRadial", "nb", "knn", "gbm")  # CHANGE HERE!
for(m in modelos) {
  for(i in seq_along(indices2)) {
    cols_df <- data.frame()
    models_boot[[paste0("boot",i)]] <- list()
    for(c in 1:nrow(cols)) {
      #boot_result <- tryCatch(calculate_accuracy(data=B4_top_genes, indices=indices2[[i]], cols=as.numeric(cols[c,]), modelo=m),error = function(e) NA)
      boot_result <- calculate_accuracy(data=B4_top_genes, indices=indices2[[i]], cols=as.numeric(cols[c,]), modelo=m)
      if(length(boot_result)<1) next
      models_boot[[paste0("boot",i)]][[c]] <- boot_result$model
      cols_df <- rbind(cols_df,
                       cbind(model=m,boot=i,cols=c,acc_treino=boot_result$acc_treino,
                             acc_teste=boot_result$acc_teste))
    }
    cols_df$mean_acc_test <- median(cols_df$acc_teste)
    accuracy_df <- rbind(accuracy_df,cols_df)
  }
  write.csv(accuracy_df,paste0("model_", m,".csv"))
}
