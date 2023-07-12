setwd("/home/keerthana2/AML/SL interactions/")
library(data.table)
library(dplyr)
library(stringr)
library(tidyr)

#depmap
depmap<- fread("final_depmap.csv")
depmap<- depmap[,-12]

#cgc data
cgc_genes<- fread("/home/keerthana/Downloads/Census_allMon May  8 06_07_52 2023.csv")

#filtering cgc genes in depmap
depmap<- filter(depmap, depmap$Gene.name%in%cgc_genes$`Gene Symbol`)

#read final_cgc.csv data
cgc<- fread("final_cgc.csv")
#length(unique(cgc$`Sample name`))

# no.of cell lines that are unique
cell_lines<- unique(cgc$`Sample name`)


# Step 1: Create an empty dataframe
gene_pairs <- data.frame(Gene.cgc = character(),
                         Gene.dep = character(),
                         dep.score= numeric(),
                         Cell_line = character(),
                         Mutation_type = character(),
                         Mutation_AA = character())
#loop to create SL pairs
for (i in 1:10){
  #cell_lines will go from 1 to 10
  cgc_data <- filter(cgc, cgc$`Sample name`== cell_lines[i])
  
  #using cgc dataframe to filter for cell line[1]
  filtered_gene_cgc<-cgc$`Gene name`[cgc$`Sample name` == cell_lines[i]]
  mut<-cgc$`Mutation Description`[cgc$`Sample name` == cell_lines[i]]
  mut_aa<- cgc$`Mutation AA`[cgc$`Sample name` == cell_lines[i]]
  #using depmap dataframe to filter for cell line[1]
  filtered_dep<-select(depmap,cell_lines[i]) 
  colnames(filtered_dep)<-"dep.score"
  # filtered_gene_cgc[1] has length of 1 to 51 (it'll vary according to gene-mutations)
#loop should be for length of filtered_gene_cgc
for (j in 1:length(filtered_gene_cgc)){
  pairs <- data.frame(Gene.cgc = filtered_gene_cgc[j],
                      Gene.dep = depmap$Gene.name,
                      dep.score=filtered_dep[,1],
                      Cell_line=cell_lines[i],
                      Mutation_type = mut[j],
                      Mutation_AA = mut_aa[j])
  gene_pairs <- rbind(gene_pairs, pairs)
}
}

write.table(gene_pairs,"gene_pairs.csv", col.names = T,row.names = F)

#gene_pairs_dep_sorted<-filter(gene_pairs,gene_pairs$dep.score>=0.7)
#write.table(gene_pairs_dep_sorted,"gene_pairs_sorted_dep.csv", col.names = T,row.names = F)
