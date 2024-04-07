setwd("~/AML/data")
library(oligo)
library(affy)
library(GEOquery)
library(dplyr)
library(tidyverse)

#gse matrix download
gse <- getGEO("GSE2191", GSEMatrix = TRUE)

if (length(gse) > 1) idx <- grep("GPL8300", attr(gse, "names")) else idx <- 1
gse <- gse[[idx]]
ex <- exprs(gse)

# log2 transform
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
ex <- log2(ex) }

#extract exprsn data
normalized_exp<-as.data.frame(ex)

#read gse
gse <- getGEO("GSE2191", GSEMatrix = TRUE)

# fetch feature data to get ID - gene symbol mapping
feature.data <- gse$`GSE2191_series_matrix.txt.gz`@featureData@data

# subset
###############################
feature..data<- feature.data[,c(1,7,11,12)]
###############################

normalized.expr_final <- normalized_exp%>%
  rownames_to_column(var = 'ID') %>%
  inner_join(., feature..data, by = 'ID')

#for 21 samples:::: 1,23:25,(2:22)
normalized.expr_final <- normalized.expr_final[,c(1,60:62,2:59)]

normalized.expr_final$`Gene Symbol`<- ifelse(normalized.expr_final$`Gene Symbol`=="",NA, normalized.expr_final$`Gene Symbol`)
normalized.expr_final$ENTREZ_GENE_ID<- ifelse(normalized.expr_final$ENTREZ_GENE_ID=="",NA, normalized.expr_final$ENTREZ_GENE_ID)
normalized.expr_final<- normalized.expr_final[!(is.na(normalized.expr_final$`Gene Symbol`) & is.na(normalized.expr_final$ENTREZ_GENE_ID)),]

#no.of probes
no_of_probes<-unique(normalized.expr_final$ID)
probe<-length(no_of_probes)

#no.of genes
genes<- normalized.expr_final[,c(3,4)]
genes<- strsplit(genes$`Gene Symbol`,split= '///' )
unlisting_gene<- data.frame(matrix(unlist(genes),byrow=TRUE),stringsAsFactors=FALSE)
final_gene<- nrow(unique(unlisting_gene))

# fetch pheno data
pheno.data<- gse$`GSE2191_series_matrix.txt.gz`@phenoData@data
pheno..data<- pheno.data[,c(2,6,8)]
pheno..data<- as.data.frame(t(pheno..data))

#save phenodata
write.table(pheno..data,"~/AML/data/GSE2191/pheno.txt",sep=" ",row.names=T, col.names =F)
write.table(normalized.expr_final,"~/AML/data/GSE2191/norm.txt",sep=" ",row.names=F)

pheno <- read.table("~/AML/data/GSE2191/pheno.txt")
norm <- read.table("~/AML/data/GSE2191/norm.txt", fill=T)

GSE2191<- bind_rows(pheno, norm)
write.table(GSE2191,"~/AML/data/GSE2191/Final_Matrix_GSE2191.txt",sep=" ",row.names=F, col.names = F)


