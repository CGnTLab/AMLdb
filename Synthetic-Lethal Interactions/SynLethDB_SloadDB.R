setwd("/home/keerthana2/AML/SL interactions/")
library(data.table)
library(dplyr)
library(stringr)
library(tidyr)


#read human SL data from SynLethDB
human_sl<- fread("Human_SL.csv")
human_sl<- human_sl[,c(1,3,5:7)]
colnames(human_sl)[c(1,2)]<-c("Gene.cgc","Gene.dep")

#read gene pairs data
gene_pairs<- fread("gene_pairs.csv")

#merging both
merged<- merge(gene_pairs,human_sl,by=c("Gene.cgc","Gene.dep"))

#renaming columns
colnames(merged)[c(7:9)]<-c("SynLethDB_cell_line","PubMed_ID","SynLethDB_source")

#assigning Yes to the column values
merged$SynLethDB<- "Yes"


#non_sl_pairs<-fread("Human_nonSL.csv")
#human_non_sl<- non_sl_pairs[,c(1,3,5:7)]
#colnames(human_non_sl)[c(1,2)]<-c("Gene.cgc","Gene.dep") 
#non_sl_merged<- merge(gene_pairs,human_non_sl,by=c("Gene.cgc","Gene.dep") )

#joining both gene pairs and merged file
SL_plus_remaining<- full_join(gene_pairs, merged)

#assigning No wherever gene pair has not been reported in SynLethDB
SL_plus_remaining$SynLethDB[SL_plus_remaining$SynLethDB==""]<-NA
SL_plus_remaining$SynLethDB[is.na(SL_plus_remaining$SynLethDB)]<-"No"

#keep unique cell line names 
SL_plus_remaining$SynLethDB_cell_line <- sapply(strsplit(SL_plus_remaining$SynLethDB_cell_line , ";"), function(x) paste(unique(x), collapse = ";"))


#SLOAD 
sload_aml<- fread("LAML_RESULT.TXT")
sload_aml<-sload_aml[,-3]
colnames(sload_aml)[c(1:2)]<-c("Gene.cgc","Gene.dep")

#merge gene pairs with SLOAD data
merged2<- merge(gene_pairs,sload_aml,by=c("Gene.cgc","Gene.dep"))

#assigning Yes to the column values
merged2$SloadDB<- "Yes"

#joining both gene pairs and merged file
SL_plus_remaining_sload<- full_join(SL_plus_remaining,merged2)

#assigning No wherever gene pair has not been reported in SLOAD
SL_plus_remaining_sload$SloadDB[SL_plus_remaining_sload$SloadDB==""]<-NA
SL_plus_remaining_sload$SloadDB[is.na(SL_plus_remaining_sload$SloadDB)]<-"No"


#save the file
write.csv(SL_plus_remaining_sload,"/home/keerthana2/AML/SL interactions/SynLethDB_SloadDB.csv", row.names = F)
