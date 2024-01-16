# Read in the file. Don't set row names yet
# Note if using R < 4.0.0, set stringsAsFactors = FALSE in read.delim
data <- read.delim("C:/Users/Ryan/Downloads/gene_expected_count.annot.txt", row.names = NULL)
# Deal with genes that don't have annotated gene symbols (external_gene_name)
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
