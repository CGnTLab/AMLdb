library(data.table)
library(dplyr)
library(stringr)
library(tidyr)
setwd("/home/keerthana2/AML/SL interactions/")


#read the cgc file
cgc<- fread("/home/keerthana/Downloads/Census_allMon May  8 06_07_52 2023.csv")

cosmic<-fread("/home/keerthana/Downloads/CosmicCLP_MutantExport (1).tsv/CosmicCLP_MutantExport.tsv")

#test<- head(cosmic)
cosmic<-cosmic[,c(1,5,22,21)]

#Remove '_' from Gene name 
cosmic$`Gene name`<- sub("_.*", "", cosmic$`Gene name`)

#common genes betwenn cosmic and cgc
cgc_gene<- cgc[,1]
colnames(cgc_gene)<-"Gene name"
final_cgc<- unique(merge(cgc_gene, cosmic, by="Gene name"))
length(unique(final_cgc$`Gene name`))

#filter silent, unknown, &synonymous mutations
final_cgc2<-filter(final_cgc, !(final_cgc$`Mutation Description`=="Unknown" | final_cgc$`Mutation Description`=="Substitution - coding silent" 
                                   | final_cgc$`Mutation Description`=="Nonstop extension"))

#remove '-' from Sample name
final_cgc2$`Sample name`<-gsub("-", "", final_cgc2$`Sample name`)
length(unique(final_cgc2$`Gene name`))

#filetring the common cell lines
depmap<- fread("final_depmap.csv")
cell_lines<- colnames(depmap)[2:11]
final_cgc3<- filter(final_cgc2, final_cgc2$`Sample name`%in% cell_lines)
length(unique(final_cgc3$`Sample name`))
length(unique(final_cgc3$`Gene name`))

write.csv(final_cgc3,"final_cgc.csv",row.names = F)

