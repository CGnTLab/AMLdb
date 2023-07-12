library(data.table)
library(dplyr)
library(stringr)
library(janitor)
library(tidyr)
setwd("/home/keerthana2/AML/SL interactions/")


#read depmap data a
dep<- fread("/home/keerthana2/AML/DEPMAP/crispr_gene_dependency_final.csv")
dep<- dep[-c(2:6),]
dep<- row_to_names(dep,1)
colnames(dep)[1]<-"Gene.name"
write.csv(dep, "/home/keerthana2/AML/SL interactions/depmap.csv", row.names = F)


#read cosmic complete mutation data
cos_mut_mod<- fread("/home/keerthana2/AML/cosmic_aml.csv")
cos_mut_mod<-unique(cos_mut_mod)
length(unique(cos_mut_mod$`Gene name`))
length(unique(cos_mut_mod$`Cell line`))

cos_mut_mod$`Cell line`<- gsub("-", "", cos_mut_mod$`Cell line`)
cos_mut_mod$`Gene name`<- sub("_.*", "", cos_mut_mod$`Gene name`)
colnames(cos_mut_mod)[1]<-"Gene.name"

#filter silent, unknown, &synonymous mutations
cos_mut_mod<-filter(cos_mut_mod, !(cos_mut_mod$`Mutation Description`=="Unknown" | cos_mut_mod$`Mutation Description`=="Substitution - coding silent" 
                                    | cos_mut_mod$`Mutation Description`=="Nonstop extension"))

#length(unique(cos_mut_mod$`Gene.name`))
#length(unique(cos_mut_mod$`Cell line`))
write.csv(cos_mut_mod,"cmd_aml.csv", row.names = F)


#CMD-DepMap
c_dep<-colnames(dep)[2:27]
c_cosmic<-unique(cos_mut_mod$`Cell line`)

common_cell_lines<-c_cosmic[c_cosmic%in%c_dep]
#common_cell_lines<-unique(cos_mut_mod$`Cell line`[cos_mut_mod$`Cell line`%in%colnames(dep)])

final_depmap<-select(dep,common_cell_lines)
final_depmap$Gene.name<-dep$Gene.name
final_depmap<-final_depmap[,c(11,1:10)]

#final_depmap<-fread("final_depmap.csv")

#count the number of column values >0.7 rowise
final_depmap$counts <- rowSums(final_depmap[,c(2:11)] > 0.7)

#filter values of counts>=5
final_depmap<-filter(final_depmap,final_depmap$counts==10)
write.csv(final_depmap,"final_depmap.csv", row.names = F)
