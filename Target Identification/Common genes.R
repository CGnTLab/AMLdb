library(dplyr)
library(data.table)
# Set the path to the folder containing the CSV files
path <- "/home/keerthana/AML/GENE_ANALYSIS/QN/"

# Get a list of all CSV files in the folder
csv_files <- list.files(path, pattern = ".csv")

# Create an empty list to store the dataframes
df_list <- lapply(csv_files, function(file) {
  filepath <- file.path(path, file)
  read.csv(filepath)
})

# Use a loop to unlist each dataframe and save as a new object
for (i in 1:length(df_list)) {
  df_unlisted <- df_list[[i]]
  assign(paste0("df_", i), df_unlisted)
}

# Extract the Gene.Symbol column from each dataframe and store in a list
gene_symbol_list <- lapply(df_list, function(x) x$Gene.Symbol)

# Use Reduce() and intersect() to find common gene symbols across all dataframes
common_gene_symbols <- Reduce(intersect, gene_symbol_list)

# common_gene_symbols will contain the common gene symbols across all dataframes
write.table(common_gene_symbols,"/home/keerthana/AML/GENE_ANALYSIS/QN/common_gene_symbols.txt",sep=",",row.names = F,col.names = F,quote = F)

#crispr_commonessentialgenes
common_essential_genes<- fread("~/Downloads/CRISPR_common_essentials.csv")
common_essential_genes<- common_essential_genes[,-2]
colnames(common_essential_genes)<- "Gene.Symbol"
common_genes<- read.table("/home/keerthana/AML/GENE_ANALYSIS/QN/common_gene_symbols.txt")
colnames(common_genes)<- "Gene.Symbol"

#removing common essential genes from the gene list
aml_genes<-filter(common_genes,!(common_genes$Gene.Symbol %in% common_essential_genes$Gene.Symbol))

#save the file
write.table(aml_genes,"/home/keerthana/AML/GENE_ANALYSIS/QN/aml_genes.txt",sep=",",row.names = F,col.names = F,quote = F)
