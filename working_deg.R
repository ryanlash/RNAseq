setwd("C:/Users/rylash/OneDrive - Michigan Medicine/Desktop/RNAseq_anals")
#gene_expected<- read.table('gene_expected_count.csv', header=TRUE, sep=',', fill=TRUE)
#gene_tpm<-read.table('gene_tpm_annot.csv', header=TRUE, sep=',', fill=TRUE)
#gene_fpkm<-read.table('gene_fpkm_annot.csv', header=TRUE, sep=',', fill=TRUE)


###### THIS IS CODE FROM THE AGC FINAL REPORT FOUND IN 'report_final.html'#######
data <- read.delim("gene_expected_count.csv", sep=',', row.names = NULL)# Deal with genes that don't have annotated gene symbols (external_gene_name)
# Use ENSEMBL ID if gene symbol not available
data$external_gene_name <- ifelse(
  data$external_gene_name == ".",
  data$gene_id,
  data$external_gene_name
)
# Deal with duplicated gene symbols
# Combine gene symbol with ENSEMBL ID if non-unique
data$external_gene_name <- ifelse(
  duplicated(data$external_gene_name),
  paste(data$external_gene_name, data$gene_id, sep="_"),
  data$external_gene_name
)
# Then we can use the gene symbol column as the row names,
# and subset the count data for further analysis
rownames(data) <- data$external_gene_name
count.data <- data[,5:ncol(data)] # All columns after 4 are count data

########## Subset counts data to include only useful samples ####################
count.data2<-count.data[,c('X9821.RL.1', 'X9821.RL.3', 'X9821.RL.5', 'X9821.RL.6')]
View(count.data2)
#create column data frame with treatment information for deseq2
col_data2 <- data.frame(
  sample_ID = c('X9821.RL.1', 'X9821.RL.3', 'X9821.RL.5', 'X9821.RL.6'),
  treatment = c('sham', 'lesion', 'sham', 'lesion')
)
#treatments need to be designated as factors
col_data2$treatment <- factor(col_data2$treatment)
col_data2$treatment <- relevel(col_data2$treatment, ref = "sham")
# ########Running Sleuth
# library(sleuth)
# sample_to_covariates <- list(
#   X9821.RL.1 = c(treatment = "sham"),
#   X9821.RL.3 = c(treatment = "lesion"),
#   X9821.RL.5 = c(treatment = "sham"),
#   X9821.RL.6 = c(treatment = "lesion")
# )
# col_data2$treatment <- factor(col_data2$treatment)
# 
# # Prepare count data with sleuth
# so <- sleuth_prep(count_data = count.data2, 
#                   sample_table = col_data2, 
#                   sample_to_covariates = sample_to_covariates, 
#                   gene_mode = TRUE)
# so <- sleuth_prep(count_data = count.data2, 
#                   sample_table = col_data2, 
#                   gene_mode = TRUE)

#Running DESeq2
library(DESeq2)
dds2 <- DESeqDataSetFromMatrix(countData = count.data2,
                              colData = col_data2,
                              design = ~ treatment)

dds2<-DESeq(dds2)
results_deseq<-results(dds2)
#dds$design <- ~ treatment
BiocManager::install('apeglm')
write.csv(as.data.frame(significant_genes_p1), 
          file="results_deseq_p1.csv")
# res01 <- results(dds2, alpha=0.01)
# summary(res01)
# res001 <- results(dds2, alpha=0.001)
# summary(res001)
resLFC <- lfcShrink(dds2, coef="treatment_sham_vs_lesion", type="apeglm")
resLFC
plotMA(resLFC)
plotMA(results_deseq)
plotMA(results_deseq,  xlab='mean of normalized counts', ylim(-8,8), colNonSig='navy',colSig='orangered')
plotCounts(dds2, gene=which.min(results_deseq$padj), intgroup="condition")
################# Create Volcano Plot############################################
with(results_deseq, plot(log2FoldChange, -log10(padj),
                         main="Volcano Plot",
                         xlim=c(-10, 10),
                         ylim=c(0, 8),
                         pch=20,
                         col=ifelse(padj < 0.05, "red", "black")))
library(ggplot2)
results_df<-as.data.frame(results_deseq)
ggplot(results_df, aes(x = log2FoldChange, y = -log10(padj), color = factor(results_df$padj < 0.05))) +
  geom_point(size = 1.5, alpha = 0.7) +
  scale_color_manual(values = c("navy", "orangered")) +
  labs(x = "Log2Fold Change",
       y = "-Log10(Adjusted p-value)",
       color = "Significant (padj<0.05)") +
  theme_light()+
  xlim(-10,10)+
  ylim(0,15)+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "hotpink") +  
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "red")+
  geom_hline(yintercept = -log10(0.001), linetype = "dashed", color = "darkred")+# Add horizontal line at -log10(0.05) 
  geom_text(data = subset(results_df, rownames(results_df) %in% c("Cacna1h", 'Slc17a7', "Sst",'Gabbr2')), 
            aes(label = rownames(subset(results_df, rownames(results_df) %in% c("Cacna1h", 'Slc17a7', "Sst",'Gabbr2')))), 
            vjust = -1, hjust = 0.4, size = 4, color = "black")
ggplot(results_deseq, aes(x = log2FoldChange, y = -log10(padj), color = factor(results_deseq$padj < 0.05))) +
  geom_point(size = 1.5, alpha = 0.7) +
  scale_color_manual(values = c("navy", "orangered")) +
  labs(x = "Log2 Fold Change",
       y = "-Log10(Adjusted p-value)",
       title = "Volcano Plot of DEGs",
       color = "Significant (padj<0.05)") +
  theme_light()+
  xlim(-10,10)+
  ylim(0,15)+  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "hotpink")   # Add horizontal line at -log10(0.05) 
################Subsetting significant genes based on p-values###################
significant_genes_p1 <- subset(results_deseq, padj < .1)
significant_genes_p05b <- subset(results_deseq, padj < .05)
significant_genes_p01 <- subset(results_deseq, padj < .01)
significant_genes_p001 <- subset(results_deseq, padj < .001)
plotMA(results_deseq)
# Add labels for significant genes ##Currently doesn't work
#significant_genes_01 <- results_deseq[results_deseq$padj < 0.01, ]
#text(significant_genes_001$log2FoldChange, -log10(significant_genes_001$padj), labels=rownames(significant_genes), col="red", cex=0.7)
#head(results_deseq)

###GENERATING SUBSET DDS BASED ONLY ON GENES WITH DESIGNATED SIGNIFICANT P VALUE###
DE_genes_names <- rownames(res001)
DE_genes_names <- rownames(significant_genes_p001)

# Subset the original DESeqDataSet based on differentially expressed genes
subset_dds <- dds2[DE_genes_names, ]
results_sub <- results(subset_dds)

############################Heatmap code#########################################
library(pheatmap)
###Subsetting for log fold change
log_fold_changes <- results_sub$log2FoldChange
gene_names <- rownames(subset_dds)

# Create a data frame with the log-fold changes and gene names
heatmap_data <- data.frame(LogFC = log_fold_changes, Gene = DE_genes_names)
pheatmap(heatmap_data)
pheatmap(heatmap_data,
         clustering_distance_cols = 'euclidean',
         main = "Heatmap of DEGs (p<0.001)")


### Assuming subset_dds is your DESeqDataSet###
heatmap_matrix <- assay(subset_dds, name=)
#heatmap_head<-head(heatmap_matrix,2:20)
# Scale the counts (optional, but can be useful for visualization)
scaled_matrix <- scale(heatmap_matrix)
#heatmap_matrix_sg2<-assay(significant_genes2)
# Create the heatmap ##NEED TO FIX THE SCALING 
pheatmap(heatmap_matrix,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         clustering_distance_cols = 'euclidean',
         color = colorRampPalette(c("navy", "white", "orangered"))(100),
         scale='row',
         fontsize_row = 8,
         main = "Heatmap of DEGs (p<0.001)")

#######PCA EXPLORER##### 

###################TSNE plots####################################################
library(Rtsne)
# Transform the DESeqDataSet using rlog or vst
dds_transformed <- rlog(dds2)  # or vst(dds)

# Extract the transformed data matrix
data_matrix <- assay(dds_transformed)

# Run t-SNE
tsne_result <- Rtsne(data_matrix, dims = 2, perplexity = 10, check_duplicates = FALSE)

# Create a data frame with the t-SNE coordinates
tsne_df <- as.data.frame(tsne_result$Y)

# Plot the t-SNE plot
plot(tsne_df, col = col_data2$treatment, pch = 16, main = "t-SNE Plot")
legend("topright", legend = levels(col_data2$treatment), col = 1:2, pch = 16, title = "Treatment")

##########################UMAP Plots ############################################
library(umap)

dds_norm <- vst(dds2)
normalized_counts <- assay(dds_norm)
normalized_counts<-t(normalized_counts)
umap_result<-umap::umap(normalized_counts, n_neighbors=2)

umap.defaults






# Assuming dds is your DESeqDataSet
dds_transformed <- rlog(dds2)  # or vst(dds)
data_matrix <- assay(dds2)
data_matrix2<-t(data_matrix)
library(umap)
# Run UMAP
umap_result <- umap(data_matrix, n_neighbors = 15, n_components = , min_dist = 0.25)

# Create a UMAP plot
umap_df <- as.data.frame(umap_result$layout)
plot(umap_df, col = col_data2$treatment, pch = 16, main = "UMAP Plot")
legend("topright", legend = levels(col_data$treatment), col = 1:2, pch = 16, title = "Treatment")

######################GO and KEGG Enrichment Analysis############################

BiocManager::install("clusterProfiler")
BiocManager::install("org.Rn.eg.db")
library(clusterProfiler)
library(org.Rn.eg.db)

# Example gene list for rats
rat_gene_list <- c("Gene1", "Gene2", "Gene3", ...)

# Perform KEGG enrichment analysis for rats
kegg_enrich_rat <- enrichKEGG(gene = significant_genes_p001, organism = "rno", keyType = "org.Rn.eg", pvalueCutoff = 0.05)
print(kegg_enrich_rat)

###############RUNNING DESeq2 on entire counts data set##########################

library('DESeq2')
str(count.data)
##set up matrix for deseq dds variable, counts data is finished
#col data needs to be created as a two level factor variable 
#and ran for construction of dds variable to complete analysis
col_data <- data.frame(
  sample_ID = c('X9821.RL.1', 'X9821.RL.2', 'X9821.RL.3', 'X9821.RL.4', 'X9821.RL.5', 'X9821.RL.6'),
  treatment = c('sham', 'lesion', 'lesion', 'sham', 'sham', 'lesion')
)
col_data$treatment <- factor(col_data$treatment)
dds <- DESeqDataSetFromMatrix(countData = count.data,
                              colData = col_data,
                              design = ~ treatment)


#dds$design <- ~ treatment
dds <- DESeq(dds)
results_deseq <- results(dds)
significant_genes2 <- subset(results_deseq, padj < .5)
plotMA(results_deseq)
if (!requireNamespace("EnhancedVolcano", quietly = TRUE)) {
  install.packages("EnhancedVolcano")
}
library(EnhancedVolcano)
plot(results_deseq, main="Volcano Plot", xlim=c(-2,2), ylim=c(0,8))
with(results_deseq, plot(log2FoldChange, -log10(padj),
                         main="Volcano Plot",
                         xlim=c(-10, 10),
                         ylim=c(0, 8),
                         pch=20,
                         col=ifelse(padj < 0.7, "red", "black")))

