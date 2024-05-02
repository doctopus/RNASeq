args <- commandArgs(trailingOnly = TRUE)
wrkg<- args[1]
input_dir<-"/mnt/beegfs/training/CITIWorkshops/RNASeq/data"
fc_dir<-paste0(input_dir,"/subread/featureCounts/")

##### creating DEseq object from the feature counts matrix
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(devtools))
suppressPackageStartupMessages(library(annotate))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(AnnotationDbi))
suppressPackageStartupMessages(library(DESeq2))
### metadata
#meta<-data.frame(Sample=c("Pt10","Pt14","Pt12","Pt13","Pt8","Pt9"),Response=c("PD","PD","PD","CR","CR","CR"))
meta<-data.frame(Sample=c("Pt10","Pt12","Pt8","Pt9"),Response=c("PD","PD","CR","CR"))
print("metadata created")
### read feature counts
fc_file<-paste0(input_dir,"/subread/featureCounts/featureCounts_0")
fc<-read.delim(fc_file,row.names=1,check.names = FALSE)
#num_samples=4
### filter genes
genes_to_keep <- rowSums(fc[,7:10])>1
fc_1<-fc[genes_to_keep,]
fc_1<-fc_1[,6:10]
fc_1$EnsembleID<-gsub("\\..*$","",rownames(fc_1))
fc_1<-fc_1[order(fc_1[,2],decreasing = TRUE),]
fc_1<-fc_1[!duplicated(fc_1$EnsembleID),]

setwd(paste0(wrkg,"/outs/counts"))
sp<-org.Hs.eg.db

### TPM normalization
geneLengths <- as_vector(subset(fc_1, select = c(Length)))
rpk <- apply( subset(fc_1, select = c(-Length,-EnsembleID)), 2, 
              function(x) x/(geneLengths/1000))
#normalize by the sample size using rpk values
tpm <- as.data.frame(apply(rpk, 2, function(x) x / sum(as.numeric(x)) * 10^6))
colSums(tpm)
tpm$GeneID<-row.names(tpm)
tpm$EnsemblID<-gsub("\\..*$", "", rownames(tpm))
tpm$GeneSymbol<- mapIds(sp, keys=tpm$EnsemblID, column=c("SYMBOL"), keytype="ENSEMBL", multiVals="first")
tpm$EntrezID<- mapIds(sp, keys=tpm$EnsemblID, column=c("ENTREZID"), keytype="ENSEMBL", multiVals="first")
rownames(tpm)<-tpm$EnsemblID
write.csv(tpm, "tpm.csv")

### counts for DESeq
fc_1$Length<-NULL
## get entrz ID and Gene symbol
fc_1$GeneSymbol<-mapIds(sp,keys=fc_1$EnsembleID, column = c("SYMBOL"), keytype = c("ENSEMBL"),multiVals = "first")
setwd(paste0(wrkg,"/outs/counts"))
saveRDS(fc_1,"clean_counts_for_DESeq.rds")
rownames(fc_1)<-fc_1$EnsembleID
write.csv(fc_1, "raw_counts_matrix_clean.csv")

### create DESeq object
library(DESeq2)
comparison<-"Response"
deseq_formula<-as.formula(paste("~",comparison,colaplse=""))
#fc_2<-fc_1[,1:6]
fc_2<-fc_1[,meta$Sample]
dds<-DESeq(DESeqDataSetFromMatrix(countData=fc_2,colData = meta,design = deseq_formula))
saveRDS(dds,"DESeqObject.rds")
## normalize count data
dds_normalized<-as.data.frame(counts(dds, normalized=TRUE))
dds_normalized<-cbind(dds_normalized,fc_1[,c("EnsembleID","GeneSymbol")])
write.csv(dds_normalized,"DESeqNormalizedDataset.csv")

### save the DESeq object as well and use that for further analysis
saveRDS(dds_normalized,"DESeqNormalizedObject.rds")
