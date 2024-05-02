args <- commandArgs(trailingOnly = TRUE)
wrkg<- args[1]
human_hall_file<- args[2]

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(devtools))
suppressPackageStartupMessages(library(annotate))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(AnnotationDbi))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(EnhancedVolcano))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(GSVA))
suppressPackageStartupMessages(library(DOSE))
suppressPackageStartupMessages(library(xCell))

## change working directory 
setwd(wrkg)
outsdir<-paste0(wrkg,"/outs")

saveFigure<-function(figure, fileName,h=7,w=7){
  pdf(file=paste0(outsdir,"/",fileName,".pdf"),height=h,width=w)
  print(figure)
  dev.off()
}

deseqobj<-paste0(wrkg,"/outs/counts/DESeqObject.rds")


#meta<-data.frame(Sample=c("Pt10","Pt14","Pt12","Pt13","Pt8","Pt9"),Response=c("PD","PD","PD","CR","CR","CR"))
meta<-data.frame(Sample=c("Pt10","Pt12","Pt8","Pt9"),Response=c("PD","PD","CR","CR"))
rownames(meta)<-meta$Sample

### change these values according to metadata and comparisons required
species="Human"
ComparisonColumn<-"Response"
factor1<-"CR"
factor2<-"PD"

sp<-org.Hs.eg.db
colors <- colorRampPalette(c("blue","white","red"))(99)

## use DESeq object
dds<-readRDS(deseqobj)
rld<-vst(dds)

## creating distance matrix
sampleDists_subset <- as.matrix(dist(t(assay(rld))))
hm<-pheatmap::pheatmap(as.matrix(sampleDists_subset),annotation_col = meta, col=colors,annotation_legend=TRUE)
saveFigure(figure=hm,fileName="DistanceMatrix",h=7,w=7)

## creating PCA plot
pca<-plotPCA(rld, intgroup="Response")
pca<-pca + geom_text(aes(label=name),vjust=2, size = 3)
saveFigure(figure=pca,fileName="PCAPlot",h=7,w=7)

### plotting DEG
e<-as.character(c(ComparisonColumn, factor1,factor2))
res<-lfcShrink(dds,contrast=e,type="normal") ## or use res<-results(dds,contrast=e)

resdata_subset <- merge(as.data.frame(res), as.data.frame(assay(rld)), by="row.names", sort=FALSE)
write.csv(resdata_subset, paste0(outsdir,"/","DifferentialExpressionAnalysis", factor1, "_vs_", factor2,".csv"),row.names = FALSE)
names(resdata_subset)[1] <- "Gene"
resdata_subset <- resdata_subset[order(resdata_subset$padj),]
resdata_subset <- resdata_subset[!is.na(resdata_subset$padj),]
#resdata_subset$EnsembleID <- gsub("\\..*$", "", resdata_subset$Gene)
resdata_subset$GeneSymbol<- mapIds(sp, keys=resdata_subset$Gene, column=c("SYMBOL"), keytype="ENSEMBL", multiVals="first")
resdata_subset$EntrezID<- mapIds(sp, keys=resdata_subset$Gene, column=c("ENTREZID"), keytype="ENSEMBL", multiVals="first")
write.csv(resdata_subset, paste0(outsdir,"/","DifferentialExpressionAnalysis_cleaned", factor1, "_vs_", factor2,".csv"),row.names = FALSE)

resdatagenes_subset <- resdata_subset[complete.cases(resdata_subset$GeneSymbol),]
resdatagenes_subset <- resdatagenes_subset[!duplicated(resdatagenes_subset$GeneSymbol), ]

resdatagenes_subset <- resdatagenes_subset[order(resdatagenes_subset$log2FoldChange),]
Top50genesdown_subset<-head(resdatagenes_subset, 50)
Top50genesup_subset<-tail(resdatagenes_subset, 50)
Top_subset<-rbind(Top50genesdown_subset, Top50genesup_subset)
#subset the counts
L<-as.character(meta$Sample)
L<-c(L, "GeneSymbol")
Top_subset<-Top_subset%>% dplyr::select(all_of(L))
TopD_subset<-data.frame(Top_subset, check.names = FALSE)
#name the rows as genesymbols
rownames(TopD_subset)<-Top_subset$GeneSymbol

## top 50 fold changes
hm<-pheatmap::pheatmap(TopD_subset[1:(length(TopD_subset)-1)],scale="row", annotation_col=meta,
                       annotation_legend =TRUE, color= colorRampPalette(c("blue","white","red"))(99), fontsize_row = 8, main="top50FoldChange")
saveFigure(figure=hm,fileName="Top50FoldChange_heatmap",h=12,w=12)

## variable genes
Ds2<-resdatagenes_subset%>% dplyr::select(all_of(L), "GeneSymbol")
Ds2<-Ds2[!duplicated(Ds2$GeneSymbol), ]
rownames(Ds2)<-Ds2$GeneSymbol
Ds2<-Ds2%>% dplyr::select(-"GeneSymbol")
#top 100 variable genes 
topVarGenes <- head(order(-genefilter::rowVars(Ds2)),100)
mat <- Ds2[topVarGenes, ]
mat <- mat - rowMeans(mat)

#plot the variable genes in heatmap
hm<-pheatmap::pheatmap(mat,color= colorRampPalette(c("blue","white","red"))(99), annotation_col=meta,
                       annotation_legend =TRUE, scale="row", fontsize_row = 8, show_rownames=T, main="top100variablegenes")
saveFigure(figure=hm,fileName="Top100VariableGenes",h=12,w=12)


## volcano plot
vp<-EnhancedVolcano(resdata_subset,
                    lab = resdata_subset$GeneSymbol,x = 'log2FoldChange',
                    y = 'pvalue',
                    pCutoff = 10e-4,
                    FCcutoff = 1)
saveFigure(figure=vp,fileName="VolcanoPlot",h=8,w=8)



### pathway analysis with hallmark geneset
### create a resources directory 
hall<-read.gmt(human_hall_file)


resdatagenes_gsea <- resdatagenes_subset
resdatagenes_gsea <- resdatagenes_gsea[!is.na(resdatagenes_gsea$log2FoldChange),]
resdatagenes_gsea <- resdatagenes_gsea[order(-resdatagenes_gsea$log2FoldChange),]
resdatagenes_gsea$FC_pval <- (resdatagenes_gsea$log2FoldChange)
logFC.l2n <- resdatagenes_gsea[order(-resdatagenes_gsea$FC_pval),]$FC_pval
names(logFC.l2n) <- resdatagenes_gsea[order(-resdatagenes_gsea$FC_pval),]$GeneSymbol

gsea.hall.l2n <- GSEA(logFC.l2n, TERM2GENE=hall, verbose=FALSE, pvalueCutoff=1)
gsea.hall.l2n.df <- as.data.frame(gsea.hall.l2n@result)
write.csv(gsea.hall.l2n.df, paste0(outsdir,"/","GSEA_output", factor1, "_vs_", factor2,".csv"),row.names = FALSE)

dp<-dotplot(gsea.hall.l2n, x="NES", showCategory=50, orderBy= "NES",font.size = 7) 
saveFigure(figure=dp,fileName="HallmarkPathwayAnalysis")

## Can be run from counts. Dont need DESeq object for this. Cannot be run on cell lines
counts<-paste0(wrkg,"/outs/counts/tpm.csv")
Df<-read.csv(counts, check.names = FALSE)
Df<-Df[!duplicated(Df$GeneSymbol),]
Df<-Df[!is.na(Df$GeneSymbol),]
rownames(Df)<-Df$GeneSymbol
Df<-Df[,meta$Sample]

### xCell

xOut<-xCellAnalysis(Df)
xOut<-data.frame(xOut, check.names = FALSE)

write.csv(xOut, paste0(outsdir, "/XCellScores.csv"),row.names = TRUE)

xOut$avg<-rowSums(xOut)/ncol(xOut)
xOut<-xOut[which(xOut$avg != xOut$Pt10),]
xOut$avg<-NULL
rownames(meta)<-meta$Sample

hm<-pheatmap::pheatmap(xOut,
                   annotation = meta,
                   show_rownames=T,
                   cluster_rows = T, 
                   cluster_cols = T, 
                   scale="row",
                   color=colorRampPalette(c("blue", "white", "red"))(99))
saveFigure(figure=hm,fileName="xCell_heatmap")



