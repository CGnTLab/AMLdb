setwd("/home/keerthana2/AML/SL interactions/")
library(data.table)
library(dplyr)
library(stringr)
library(tidyr)

#SynLeth_sload
syn_slo<- fread("SynLethDB_SloadDB.csv")

#READ ccle
ccle<- fread("CCLE_expression.csv")
#colnames(ccle)<- sub("\\s*\\(.*\\)", "", colnames(ccle)[2:19222])

#need only gene names
colnames(ccle)[2:19222]<- sub("\\s*\\(.*\\)", "", colnames(ccle)[2:19222])
colnames(ccle)[1]<- "DepMap_ID"


#read depmap data 
dep<- fread("/home/keerthana2/AML/DEPMAP/crispr_gene_dependency_final.csv")
dep<- dep[1,]
depId<- colnames(dep)[-1]
celllines<-dep[1,2:27]
df<-as.data.frame(t(celllines))
df<- tibble::rownames_to_column(df, "DepMapID")
final_depmap<- fread("final_depmap.csv")
cell_lines<- colnames(final_depmap)[2:11]
df<- filter(df, df$V1%in% cell_lines)
df_10<- df[,1]
 
#map ccle DepMap IDs  with the 10 common cell lines in DepMap
ccle_final<- filter(ccle, ccle$DepMap_ID%in%df_10)
ccle_final<- as.data.frame(t(ccle_final))
ccle_final<- tibble::rownames_to_column(ccle_final,"DepMap_ID")
colnames(ccle_final)<- ccle_final[1,]
colnames(ccle_final)[2:11]<- cell_lines
ccle_final<- ccle_final[-1,]

#filtering those genes with expression>=3 
ccle_final$counts <- rowSums(ccle_final[,c(2:11)] >=3)
 
#filter values of counts>=3 in all 10 cell lines
ccle_final<-filter(ccle_final,ccle_final$counts==10)
write.csv(ccle_final,"final_ccle.csv", row.names = F) 
