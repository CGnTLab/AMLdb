library(data.table)
library(dplyr)
setwd("/home/keerthana2/AML/SL interactions/")

#read the SL data
df<- fread("SynLethDB_SloadDB.csv")

#read final_ccle.csv
ccle<- fread("final_ccle.csv")
ccle<- ccle[,-12]


###for cgc genes and cell lines

#make empty datafarme
result1<- data.frame(Gene.cgc = character(),
                    Cell_line = character(),
                    Exp.cgc = character(),
                    Exp.cgc.value = numeric())

genes.cgc<- df$Gene.cgc
cell_lines<-df$Cell_line
#loop to ad expression status and values for cgc genes
for(i in 1:18614){
  Cell_line<-cell_lines[i]
  if(genes.cgc[i]%in%ccle$DepMap_ID){
    Exp.cgc = "Yes"
    filtered_cgc<- filter(ccle, ccle$DepMap_ID==genes.cgc[i])
    exp_value_cgc<- select(filtered_cgc, cell_lines[i])
    colnames(exp_value_cgc)<- "Exp.cgc.value"
    result_cgc<- data.frame(Gene.cgc = genes.cgc[i],
                            Cell_line = cell_lines[i],
                            Exp.cgc = Exp.cgc,
                            Exp.cgc.value = exp_value_cgc)
  }else{
    Exp.cgc = "No"
    Exp.cgc.value = "NA"
    result_cgc<- data.frame(Gene.cgc = genes.cgc[i],
                            Cell_line = cell_lines[i],
                            Exp.cgc = Exp.cgc,
                            Exp.cgc.value = Exp.cgc.value)
  }
  result1<- rbind(result1,result_cgc)
}



###for dep genes and cell lines
#make empty datafarme
#make empty datafarme
result2<- data.frame(Gene.dep = character(),
                    Cell_line = character(),
                    Exp.dep = character(),
                    Exp.dep.value = numeric())

genes.dep<- df$Gene.dep
cell_lines<-df$Cell_line
#loop to ad expression status and values for depmap genes
for(i in 1:18614){
  Cell_line<-cell_lines[i]
  if(genes.dep[i]%in%ccle$DepMap_ID){
    Exp.dep = "Yes"
    filtered_cgc<- filter(ccle, ccle$DepMap_ID==genes.dep[i])
    exp_value_dep<- select(filtered_cgc, cell_lines[i])
    colnames(exp_value_dep)<- "Exp.dep.value"
    result_dep<- data.frame(Gene.dep = genes.dep[i],
                            Cell_line = cell_lines[i],
                            Exp.dep = Exp.dep,
                            Exp.dep.value = exp_value_dep)
  }else{
    Exp.dep = "No"
    Exp.dep.value = "NA"
    result_dep<- data.frame(Gene.dep = genes.dep[i],
                            Cell_line = cell_lines[i],
                            Exp.dep = Exp.dep,
                            Exp.dep.value = Exp.dep.value)
  }
  result2<- rbind(result2,result_dep)
}

##combining result1 &2
final_df<- cbind(result1,result2)
final_df<-final_df[,-6]
final_df<- final_df[,c(1,5,2,3,4,6,7)]

##recheck if the order is coorect
identical(final_df$Gene.cgc, df$Gene.cgc) & identical(final_df$Gene.dep, df$Gene.dep) & identical(final_df$Cell_line,df$Cell_line)


##combine other columns
final_df<- cbind(final_df,df[,c(3,5:11)])
identical(final_df$dep.score,df$dep.score)
identical(final_df$Mutation_type,df$Mutation_type)
identical(final_df$Mutation_AA,df$Mutation_AA)
complete_SLI<-final_df[,c(1,4,5,2,6:8,3,9,13,14,15,16:18)]
complete_SLI<-complete_SLI[,c(1,4,7,8,10,9,14,15,11,12,13,2,3,5,6)]
##save the file
write.csv(complete_SLI,"complete_SLI.csv",row.names = F)

###genes with SL
sli<- fread("complete_SLI.csv")
sl_genes<- sli[,c(1,2)]
genes<-unique(sl_genes$Gene.dep)
write.table(genes,"interaction_network_genes.txt", row.names = F, col.names = F)
