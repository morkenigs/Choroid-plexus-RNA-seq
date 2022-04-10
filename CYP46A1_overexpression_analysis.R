#Library Loading
library("BiocParallel")
library(DESeq2)
library(FactoMineR)
library(pheatmap) 
library(apeglm)
library(ggplot2)
library(Hmisc)
library(stringr)
library(ggpubr)

#Accessory functions
string.to.colors = function (string, colors = NULL) 
{
  if (is.factor(string)) {
    string = as.character(string)
  }
  if (!is.null(colors)) {
    if (length(colors) != length(unique(string))) {
      (break)("The number of colors must be equal to the number of unique elements.")
    }
    else {
      conv = cbind(unique(string), colors)
    }
  }
  else {
    conv = cbind(unique(string), rainbow(length(unique(string))))
  }
  unlist(lapply(string, FUN = function(x) {
    conv[which(conv[, 1] == x), 2]
  }))
}

setwd("/home/labs/amit/mork/R/Experiments/100121_invivo_cyp46a1")

##I) Loading/annotation/QC of the  data - don't forget to delete columns before gene annotation.
count_table=read.delim("20220103_exons_raw.txt",header = T,row.names = 1)
coldata=read.delim("samples.txt",header=T)
head(count_table)
dim(count_table)

#convert gene strings to names
gene_names=rownames(count_table)
x = strsplit(gene_names, "\\|")
x = unlist(lapply(x,FUN = function(x) {x[1]}))
rownames(count_table)=x
colnames(count_table) = as.character(t(coldata$Condition))

#Annotation
annotation_table=data.frame(coldata)
#colnames(annotation_table) = "Condition"

#QC library size
barplot(log10(colSums(count_table))) ## -> deep sequencing + homogenous data

##II)Analysis using DESeq
dds <- DESeqDataSetFromMatrix(countData = count_table,
                              colData = annotation_table,
                              design = ~ Condition)
hist(log10(1+rowSums(count_table)),n=100)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
#define reference level
dds$Condition <- relevel(dds$Condition, ref = "GFP")
#differential expression analysis
dds <- DESeq(dds)
plotDispEsts(dds)
plotMA(dds, ylim=c(-3,3),alpha = 0.01)
res <- results(dds)
resultsNames(dds)
#LFC shrinkage
resLFC <- lfcShrink(dds, coef="Condition_Cyp46a1_vs_GFP", type="apeglm")
plotMA(resLFC, ylim=c(-6,6),alpha = 0.01)
resLFC[order(resLFC$log2FoldChange, decreasing = T),]
table_LFC=as.matrix(resLFC)

#check that the highest-expressed genes are cell-type specific (cell type purity):
resLFC[order(resLFC$baseMean, decreasing = T),]

#log transformation - for visualization (can be done with rlog or VST)
##rlog
transformed_expression <- rlog(dds, blind=FALSE)
# assay function extracts the matrix of normalized values
transformed_expression=assay(transformed_expression) 

#VST transformation
vsd <- varianceStabilizingTransformation(dds,blind = T)
Normalized_table = assay(vsd)
Normalized_variance = apply(Normalized_table,MARGIN = 1,FUN = var)
Normalized_variance =Normalized_variance[order(Normalized_variance,decreasing = T)]
Top_2000_genes = names(Normalized_variance[1:2000])
Top_500_genes = names(Normalized_variance[1:500])

#sample to sample distance
sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$Condition, vsd$experiment, sep="-")
colnames(sampleDistMatrix) <- rownames(sampleDistMatrix)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors,
          show_rownames = T)
#sample-sample correlation:
sampleCor=cor(assay(vsd),method = "spearman")
rownames(sampleCor) <- paste(vsd$Condition, vsd$experiment, sep="-")
colnames(sampleCor) <- rownames(sampleDistMatrix)
pheatmap(sampleCor)
#plot PCA
plotPCA(vsd, intgroup=c("Condition", "mouse_number"))
# PCA analysis
PCA_analysis = PCA(t(Normalized_table[Top_2000_genes,]),scale.unit = F,graph = F)
par(las=1,bty="l")
#plot explained variance
barplot(PCA_analysis$eig[,2],names.arg = "",ylim=c(0,70),
        ylab="Variance explained (%)",xlab='PCs',cex.lab=1.3)


#export list of all genes+LFC+pval
write.table(table_LFC, "results.txt", sep="\t")

#prepare genes for bargraph
for_bar_all=as.data.frame(resLFC[which(resLFC$padj<0.05),2:3])
for_bar_all=for_bar_all[1:6,]

#bargraph
ggplot(data=for_bar_all, aes(x=reorder(rownames(for_bar_all),-log2FoldChange), y=log2FoldChange))+geom_bar(stat="identity", fill="#F09837",color="black",size=0.46,width=0.6)+
  geom_errorbar(aes(ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE), width=.2,position=position_dodge(.9),color="black",size=0.46)+ 
  ylab(expression("log"[2]*" f.c."))+
  theme(axis.text.x=element_text(size=10,  family="sans",face = "italic",color="black",angle = 45,vjust=1,hjust=1),axis.text.y=element_text(size=10,  family="sans",color="black"), axis.title =element_text(size=10,  family="sans", color="black"),axis.line = element_line(colour = 'black', size = 0.46),
        axis.ticks.length=unit(.15, "cm"),axis.ticks =element_line(color="black",size=0.46),axis.title.x = element_blank())+
  scale_y_continuous(expand = c(0,0))+geom_hline(yintercept=0, linetype="solid", color="black", size=0.46)+ylim(-3,2)
ggsave(filename = "barplot_cyp46a1.pdf", height=55, unit="mm",width=80)

#Wpre TPM:
wpre_raw=c(7.000,	10.000,	2.000,	12.000,	14.000,	13.000)
names(wpre_raw)=colnames(count_table)
transcripts=10^(log10(colSums(count_table)))
wpre_tpm=(wpre_raw/transcripts)*(10^6)

#plot library size
lib_size=as.data.frame(log10(colSums(count_table)))
lib_size$AB=annotation_table$AB
colnames(lib_size)=c("size","AB")
lib_size=merge(lib_size,annotation_table,by.x="AB", by.y="AB")

ggplot(data=lib_size, aes(x=AB, y=size, fill=Condition))+geom_bar(stat="identity",color="black",size=0.46)+xlab("Sample")+ylab(expression("log"[10]*" total UMI"))+
  theme(legend.position="none",axis.text.x=element_text(size=8,  family="sans",face = "italic",color="black",angle = 45,vjust=1,hjust=1),axis.text.y=element_text(size=8,  family="sans",color="black"), axis.title =element_text(size=8,  family="sans", color="black"),axis.line = element_line(colour = 'black', size = 0.126),
        axis.ticks.length=unit(.15, "cm"),axis.ticks =element_line(color="black",size=0.126),axis.title.x = element_blank(),strip.background = element_blank(),strip.placement = "outside", strip.text.x = element_text(size=8,  family="sans",color="black"))+scale_fill_manual(values=c("#F09837","#A5A5A5"))
ggsave(filename = "libSize.pdf", height=55, unit="mm",width=80)

#plot PCA
pdf("PCA.pdf", width=4, height=2)
plotPCA(vsd, intgroup=c( "Condition"))
PCA_plot=plotPCA(vsd, intgroup=c( "Condition"),returnData=T)
levels(PCA_plot$Condition)=c("GFP","Cyp46a1")

ggplot(data=PCA_plot, aes(x=PC1, y=PC2))+geom_point(aes(fill=Condition), shape=21, size=2)+scale_fill_manual(values=c("#A5A5A5","#F09837"))+
  xlab("PC1: 30% variance")+ylab("PC2: 21% variance")+
  theme( axis.text.x=element_text(size=10,  family="sans",color="black",vjust=1,hjust=1),axis.text.y=element_text(size=10,  family="sans",color="black"), axis.title =element_text(size=10,  family="sans", color="black"),axis.line = element_line(colour = 'black', size = 0.46),
         axis.ticks.length=unit(.15, "cm"),axis.ticks =element_line(color="black",size=0.46))
ggsave(filename = "PCA.pdf", height=55, unit="mm",width=80)
