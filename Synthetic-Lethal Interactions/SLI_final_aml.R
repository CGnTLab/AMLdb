library(data.table)
library(dplyr)

setwd("/home/keerthana2/AML/SL interactions")
#complete sl
sli<-fread("E:/M.Tech Thesis/SL interactions/complete_SLI.csv")

#remove self-interactions
sli <- subset(sli, sli$Gene.cgc != sli$Gene.dep)

#removing unwanted mutations
sli<-filter(sli,!(sli$Mutation_type=="Substitution - Missense" | sli$Mutation_type=="Insertion - In frame" |
                    sli$Mutation_type=="Deletion - In frame"))

#unique gene pairs
unique<-unique(sli[,c(1,2)])

#filter for dep=yes
dep_yes<- filter(sli,sli$Exp.dep=="Yes")

#load cgc data
#cosmic_cgc genes for AML
cosmic_aml<- fread("C:\\Users\\keert\\Desktop\\Census_allThu Jun  1 05_52_43 2023.csv")
cosmic_aml<- filter(cosmic_aml, !(cosmic_aml$`Tumour Types(Somatic)`=="salivary gland mucoepidermoid"))
cosmic_aml_genes<- cosmic_aml$`Gene Symbol`

#add this column to dep_yes, if present in gene.cgc yes else no
dep_yes$Gene.cgc.aml<-"No"
dep_yes$Gene.cgc.aml<- ifelse(dep_yes$Gene.cgc%in%cosmic_aml_genes, "Yes", dep_yes$Gene.cgc.aml)
#add this column to dep_yes, if present in gene.cgc yes else no
dep_yes$Gene.dep.aml<-"No"
dep_yes$Gene.dep.aml<- ifelse(dep_yes$Gene.dep%in%cosmic_aml_genes, "Yes", dep_yes$Gene.dep.aml)

#rearrage
dep_yes<- dep_yes[,c(1,16,2,17,3:15)]

#select oncogenes from cgc
cgc<- fread("C:\\Users\\keert\\Desktop\\Census_allThu Jun  1 05_52_43 2023.csv")
cgc<-cgc[,c(1,15)]

#oncogenes<- filter(cgc, cgc$`Role in Cancer`=="oncogene" |cgc$`Role in Cancer`=="oncogene, fusion")
#oncogenes<- oncogenes$`Gene Symbol`

#add this column to dep_yes if oncogene yes else no
#dep_yes$Role_in_cancer<- "Yes"
#dep_yes$Role_in_cancer<- ifelse(dep_yes$Gene.cgc%in%oncogenes, "Oncogene","TSG")

#select tsg from cgc
#tsg<- filter(cgc, cgc$`Role in Cancer`=="TSG" |cgc$`Role in Cancer`=="TSG, fusion" | cgc$`Role in Cancer`=="TSG, fusion" )
#oncogenes<- oncogenes$`Gene Symbol`

#unique(dep_yes$Mutation_type)

#removing unwanted mutations
#dep_yes_final<-filter(dep_yes,!(dep_yes$Mutation_type=="Substitution - Missense" | dep_yes$Mutation_type=="Insertion - In frame" |
#dep_yes$Mutation_type=="Deletion - In frame"))

#unique(dep_yes_final$Gene.cgc[dep_yes_final$Gene.cgc%in%cosmic_aml_genes])
colnames(cgc)[1]<-"Gene.cgc"
merged<-merge(dep_yes,cgc,by='Gene.cgc')

final_aml<-filter(merged,merged$Gene.cgc.aml=='Yes')

#unique gene pairs
unique<-unique(final_aml[,c(1,3)])

write.csv(final_aml,"C:\\Users\\keert\\Desktop\\SLI_final_aml.csv",row.names = F)

