#Library Loading
library("BiocParallel")
library(DESeq2)
library(FactoMineR)
library(apeglm)
library(pheatmap)
library(ggplot2)
library(Hmisc)
library(stringr)
library(ggpubr)
library("RColorBrewer")
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



setwd("/home/labs/amit/mork/R/Experiments/020221_CP_24OH_24h/integrated analysis with previous exp/")

exp20221=read.delim("20210201_exons_raw_24OH.txt",header = T,row.names = 1)
gene_names=rownames(exp20221)
x = strsplit(gene_names, "\\|")
x = unlist(lapply(x,FUN = function(x) {x[1]}))
rownames(exp20221)=x

exp140920=read.delim("20200903_exons_raw.txt",header = T,row.names = 1)
gene_names=rownames(exp140920)
x = strsplit(gene_names, "\\|")
x = unlist(lapply(x,FUN = function(x) {x[1]}))
rownames(exp140920)=x

exp141020=read.delim("20201012_exons_raw.txt",header = T,row.names = 1)
gene_names=rownames(exp141020)
x = strsplit(gene_names, "\\|")
x = unlist(lapply(x,FUN = function(x) {x[1]}))
rownames(exp141020)=x

pilot=read.delim("20200615A_exons_raw.txt",header = T,row.names = 1)
gene_names=rownames(pilot)
x = strsplit(gene_names, "\\|")
x = unlist(lapply(x,FUN = function(x) {x[1]}))
rownames(pilot)=x

count_table=merge(exp20221,exp140920, by=0)
rownames(count_table)=count_table$Row.names
count_table=count_table[,-1]
count_table=merge(count_table,exp141020, by=0)
rownames(count_table)=count_table$Row.names
count_table=count_table[,-1]
count_table=merge(count_table,pilot, by=0)
rownames(count_table)=count_table$Row.names
count_table=count_table[,-1]



##I) Loading/annotation/QC of the  data - don't forget to delete columns before gene annotation.
coldata=read.delim("samples_all.txt",header=T)
head(count_table)
dim(count_table)
colnames(count_table) = as.character(t(coldata$AB))

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
dds$Condition <- relevel(dds$Condition, ref = "DMSO")
#differential expression analysis
dds <- DESeq(dds)
plotDispEsts(dds)
plotMA(dds, ylim=c(-3,3),alpha = 0.01)
res <- results(dds)
resultsNames(dds)
#LFC shrinkage
resLFC <- lfcShrink(dds, coef="Condition_24_OH_vs_DMSO", type="apeglm")
plotMA(resLFC, ylim=c(-6,6),alpha = 0.01)
resLFC[order(resLFC$log2FoldChange, decreasing = T),]
table_LFC=as.matrix(resLFC)

#check that the highest-expressed genes are cell-type specific (culture purity):
resLFC[order(resLFC$baseMean, decreasing = T),]

#log transformation - for visualization (can be done with rlog or VST)
##rlog
transformed_expression <- rlog(dds, blind=FALSE)
# assay function extracts the matrix of normalized values
transformed_expression=assay(transformed_expression) 

gene_plot_function=function(gene) {
  par(las=1)
  u=data.frame(Gene=transformed_expression[gene,],
               Condition=factor(annotation_table$Condition,levels =as.character(unique(annotation_table$Condition))))
  plot(Gene~Condition,u,pch=21,bg="orange",cex=2,main=gene,col="black",cex.main=1.5,ylab="Gene expression (Log2)")
}

gene_plot_function("Apoe")
gene_plot_function("Abca1")
gene_plot_function("Cd36")

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
plotPCA(vsd, intgroup=c("Condition"))
# PCA analysis
PCA_analysis = PCA(t(Normalized_table[Top_2000_genes,]),scale.unit = F,graph = F)
par(las=1,bty="l")
#plot explained variance
barplot(PCA_analysis$eig[,2],names.arg = "",ylim=c(0,70),
        ylab="Variance explained (%)",xlab='PCs',cex.lab=1.3)

#check which DMSO is clustered with the 24OH
summary(PCA_analysis)
#AB525,241,242,243 (pilot),422--remove them 
new_count_table=count_table[ , -which(names(count_table) %in% c("AB525","AB422","AB241","AB242","AB243"))]
new_annotation_table=as.data.frame(annotation_table[-which(annotation_table$AB %in% c("AB525","AB422","AB241","AB242","AB243")),])


##II)Analysis using DESeq
new_dds <- DESeqDataSetFromMatrix(countData = new_count_table,
                              colData = new_annotation_table,
                              design = ~ Condition)
hist(log10(1+rowSums(new_count_table)),n=100)
keep <- rowSums(counts(new_dds)) >= 10
new_dds <- new_dds[keep,]

#plot library size
lib_size=as.data.frame(log10(colSums(new_count_table)))
lib_size$AB=rownames(lib_size)
colnames(lib_size)=c("size","AB")
lib_size=merge(lib_size,new_annotation_table,by.x="AB", by.y="AB")

ggplot(data=lib_size, aes(x=AB, y=size, fill=Condition))+geom_bar(stat="identity",color="black",size=0.46)+xlab("Sample")+ylab(expression("log"[10]*" total UMI"))+
  theme(legend.position="none",axis.text.x=element_text(size=8,  family="sans",face = "italic",color="black",angle = 45,vjust=1,hjust=1),axis.text.y=element_text(size=8,  family="sans",color="black"), axis.title =element_text(size=8,  family="sans", color="black"),axis.line = element_line(colour = 'black', size = 0.126),
        axis.ticks.length=unit(.15, "cm"),axis.ticks =element_line(color="black",size=0.126),axis.title.x = element_blank(),strip.background = element_blank(),strip.placement = "outside", strip.text.x = element_text(size=8,  family="sans",color="black"))+scale_fill_manual(values=c("#F09837","#A5A5A5"))
ggsave(filename = "libSize.pdf", height=55, unit="mm",width=80)

#define reference level
new_dds$Condition <- relevel(new_dds$Condition, ref = "DMSO")
#differential expression analysis
new_dds <- DESeq(new_dds)
plotDispEsts(new_dds)
plotMA(new_dds, ylim=c(-3,3),alpha = 0.01)
new_res <- results(new_dds)
resultsNames(new_dds)
#LFC shrinkage
new_resLFC <- lfcShrink(new_dds, coef="Condition_24_OH_vs_DMSO", type="apeglm")
plotMA(new_resLFC, ylim=c(-6,6),alpha = 0.01)
new_resLFC[order(new_resLFC$log2FoldChange, decreasing = T),]
new_table_LFC=as.matrix(new_resLFC)

#check that the highest-expressed genes are cell-type specific (culture purity):
new_resLFC[order(new_resLFC$baseMean, decreasing = T),]

#VST transformation
new_vsd <- varianceStabilizingTransformation(new_dds,blind = T)
new_Normalized_table = assay(new_vsd)
new_Normalized_variance = apply(new_Normalized_table,MARGIN = 1,FUN = var)
new_Normalized_variance =new_Normalized_variance[order(new_Normalized_variance,decreasing = T)]
new_Top_2000_genes = names(new_Normalized_variance[1:2000])
new_Top_500_genes = names(new_Normalized_variance[1:500])

write.table(new_Normalized_table,file="24OH_data_vst_normalized.txt")

#sample to sample distance
sampleCor=cor(assay(new_vsd),method = "spearman")
rownames(sampleCor) <- paste(new_vsd$Condition, new_vsd$experiment, sep="-")
colnames(sampleCor) <- rownames(sampleDistMatrix)
pheatmap(sampleCor)

sampleDists <- dist(t(assay(new_vsd)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(new_vsd$Condition, new_vsd$experiment, sep="-")
colnames(sampleDistMatrix) <- rownames(sampleDistMatrix)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors,
         show_rownames = T)
#plot PCA
pdf("PCA.pdf", width=4, height=2)
plotPCA(new_vsd, intgroup=c( "Condition"))
PCA_plot=plotPCA(new_vsd, intgroup=c( "Condition"),returnData=T)
levels(PCA_plot$Condition)=c("DMSO","24-OH")

ggplot(data=PCA_plot, aes(x=PC1, y=PC2))+geom_point(aes(fill=Condition), shape=21, size=2)+scale_fill_manual(values=c("#A5A5A5","#F09837"))+
xlab("PC1: 83% variance")+ylab("PC2: 6% variance")+
theme( axis.text.x=element_text(size=10,  family="sans",color="black",vjust=1,hjust=1),axis.text.y=element_text(size=10,  family="sans",color="black"), axis.title =element_text(size=10,  family="sans", color="black"),axis.line = element_line(colour = 'black', size = 0.46),
axis.ticks.length=unit(.15, "cm"),axis.ticks =element_line(color="black",size=0.46))
ggsave(filename = "PCA.pdf", height=55, unit="mm",width=80)

# PCA analysis
new_PCA_analysis = PCA(t(new_Normalized_table[new_Top_2000_genes,]),scale.unit = F,graph = F)
par(las=1,bty="l")
new_resLFC[order(new_resLFC$log2FoldChange, decreasing = T),]

genes_cholesterol=c("Abca1","Abcg1","Acat2","Apoe","Cyp51","Dhcr24","Dhcr7","Fdft1","Fdps",
"Hmgcr","Hmgcs1","Hsd17b7","Ldlr","Lss","Msmo1","Mvd","Nsdhl","Sc5d","Scarb1","Sqle",
"Stard4","Stard5")
genes_vesicular=c("Ap1s2","Cav3","Clasp2","Cltb","Dync1i1","Flot1","Itsn1","Myo1f")
genes_typeI=c("Ifih1","Ifit1","Ifit2","Ifit3","Ifitm2","Ifitm3","Ifnar1","Ifnar2","Ikbke","Irf1","Irf7","Irf9","Stat1","Stat2","B2m")
genes_igf=c("Igf1","Igf1r","Igfals","Igfbp3","Igfbp4","Igfbp6","Igfbp7")
genes_down_proteomics_overlap=c("B2m","Aebp1","Spp1","Igfbp3","Nbl1","Tgfbr3","Ldlr","Glrx")
genes_all=c(genes_cholesterol,genes_vesicular,genes_typeI,genes_igf,genes_down_proteomics_overlap)
for_bar_all=as.data.frame(new_resLFC[genes_all,2:3])
for_bar_all[c("Abca1","Abcg1","Apoe"),3]="Cholesterol export"
for_bar_all[c("Ldlr","Stard4","Stard5","Scarb1"),3]="Cholesterol import"
for_bar_all[c("Acat2","Cyp51","Dhcr24","Dhcr7","Fdft1","Fdps",
            "Hmgcr","Hmgcs1","Hsd17b7","Lss","Msmo1","Mvd","Nsdhl","Sc5d","Sqle"),3]="Cholesterol synthesis"
names(for_bar_all)[3]="subgroup"

#barplot vesicular
ggplot(data=for_bar_all[genes_vesicular,], aes(x=reorder(genes_vesicular,-log2FoldChange), y=log2FoldChange))+geom_bar(stat="identity", fill="#F09837",color="black",size=0.46,width=0.6)+
  geom_errorbar(aes(ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE), width=.2,position=position_dodge(.9),color="black",size=0.46)+ 
  ylab(expression("log"[2]*" f.c."))+
  theme(axis.text.x=element_text(size=10,  family="sans",face = "italic",color="black",angle = 45,vjust=1,hjust=1),axis.text.y=element_text(size=10,  family="sans",color="black"), axis.title =element_text(size=10,  family="sans", color="black"),axis.line = element_line(colour = 'black', size = 0.46),
        axis.ticks.length=unit(.15, "cm"),axis.ticks =element_line(color="black",size=0.46),axis.title.x = element_blank())+
  scale_y_continuous(expand = c(0,0)) 
ggsave(filename = "barplot_vesicular.pdf", height=55, unit="mm",width=80)

#barplot typeI
ggplot(data=for_bar_all[genes_typeI,], aes(x=reorder(genes_typeI,-log2FoldChange), y=log2FoldChange))+geom_bar(stat="identity", fill="#F09837",color="black",size=0.46,width=0.6)+
  geom_errorbar(aes(ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE), width=.2,position=position_dodge(.9),color="black",size=0.46)+ 
  ylab(expression("log"[2]*" f.c."))+
  theme(axis.text.x=element_text(size=10,  family="sans",face = "italic",color="black",angle = 45,vjust=1,hjust=1),axis.text.y=element_text(size=10,  family="sans",color="black"), axis.title =element_text(size=10,  family="sans", color="black"),axis.line = element_line(colour = 'black', size = 0.46),
        axis.ticks.length=unit(.15, "cm"),axis.ticks =element_line(color="black",size=0.46),axis.title.x = element_blank())+
  scale_y_continuous(expand = c(0.05,-0.05)) +
  geom_hline(yintercept=0, linetype="solid", color="black", size=0.46)
ggsave(filename = "barplot_typeI.pdf", height=55, unit="mm",width=80)

#barplot Igf
ggplot(data=for_bar_all[genes_igf,], aes(x=reorder(genes_igf,-log2FoldChange), y=log2FoldChange))+geom_bar(stat="identity", fill="#F09837",color="black",size=0.46,width=0.6)+
  geom_errorbar(aes(ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE), width=.2,position=position_dodge(.9),color="black",size=0.46)+ 
  ylab(expression("log"[2]*" f.c."))+
  theme(axis.text.x=element_text(size=10,  family="sans",face = "italic",color="black",angle = 45,vjust=1,hjust=1),axis.text.y=element_text(size=10,  family="sans",color="black"), axis.title =element_text(size=10,  family="sans", color="black"),axis.line = element_line(colour = 'black', size = 0.46),
        axis.ticks.length=unit(.15, "cm"),axis.ticks =element_line(color="black",size=0.46),axis.title.x = element_blank())+
  scale_y_continuous(expand = c(0.05,-0.05)) +
  geom_hline(yintercept=0, linetype="solid", color="black", size=0.46)
ggsave(filename = "barplot_igf.pdf", height=55, unit="mm",width=80)

#barplot cholesterol
ggplot(data=for_bar_all[genes_cholesterol,], aes(x=reorder(genes_cholesterol,-log2FoldChange), y=log2FoldChange))+geom_bar(stat="identity", fill="#F09837",color="black",size=0.46,width=0.6)+
  geom_errorbar(aes(ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE), width=.2,position=position_dodge(.9),color="black",size=0.46)+ 
  ylab(expression("log"[2]*" f.c."))+
  theme(axis.text.x=element_text(size=10,  family="sans",face = "italic",color="black",angle = 45,vjust=1,hjust=1),axis.text.y=element_text(size=10,  family="sans",color="black"), axis.title =element_text(size=10,  family="sans", color="black"),axis.line = element_line(colour = 'black', size = 0.46),
        axis.ticks.length=unit(.15, "cm"),axis.ticks =element_line(color="black",size=0.46),axis.title.x = element_blank(),strip.background = element_blank(),strip.placement = "outside", strip.text.x = element_text(size=10,  family="sans",color="black"))+
  scale_y_continuous(expand = c(0,0.3)) +facet_grid(~subgroup, scales = "free_x",space="free_x",switch="x",labeller = labeller(subgroup = label_wrap_gen(10))) +geom_hline(yintercept=0, linetype="solid", color="black", size=0.46)
ggsave(filename = "barplot_cholesterol.pdf", height=80, unit="mm",width=150)

#Prepare data for plotting GO terms from metascape results sheet
eval_vec=Vectorize(eval, vectorize.args="expr")
GO_down_24=read.delim("metascape_result_down.txt",header = T,stringsAsFactors = F)
GO_down_24=GO_down_24[grep("*Summary*",GO_down_24$GroupID),]
GO_up_24=read.delim(file = "metascape_result_up.txt",header = T,stringsAsFactors = F)
GO_up_24=GO_up_24[grep("*Summary*",GO_up_24$GroupID),]

GO_down_24$LogP_pos=abs(GO_down_24$LogP)
GO_up_24$LogP_pos=abs(GO_up_24$LogP)
GO_down_24$ratio=eval_vec(parse(text = GO_down_24$InTerm_InList))
GO_up_24$ratio=eval_vec(parse(text = GO_up_24$InTerm_InList))
GO_up_24$counts=as.numeric(lapply(GO_up_24$InTerm_InList,function(x) {unlist(strsplit(x,"/"))[1]}))
GO_down_24$counts=as.numeric(lapply(GO_down_24$InTerm_InList,function(x) {unlist(strsplit(x,"/"))[1]}))
GO_down_24$Description=capitalize(GO_down_24$Description)
GO_up_24$Description=capitalize(GO_up_24$Description)
GO_down_24$Description[11]=str_to_sentence(GO_down_24$Description[11])
GO_up_24$Description[1]=str_to_sentence(GO_up_24$Description[1])
GO_up_24$Description[13]="Homologous DNA pairing and strand exchange"

#Bubble plot of GO terms 
ggplot(data=GO_down_24[1:15,], aes(x=counts, y=reorder(Description,counts),size=ratio,color=LogP_pos))+geom_point()+
  scale_size(range = c(0.5, 3.5), name="Gene Ratio",breaks = c(0.2,0.3,0.4))+
  theme(axis.text.x=element_text(size=6.5,  family="sans",color="black"),axis.text.y=element_text(size=6.5,  family="sans",color="black"), 
        axis.title =element_text(size=6.5,  family="sans", color="black"),axis.line = element_line(colour = 'black', size = 0.46),axis.ticks.length=unit(.15, "cm"),
        axis.ticks =element_line(color="black",size=0.46),axis.title.y=element_blank(),panel.grid.major=element_line(size=0.46),legend.title = element_text(size=6.5),legend.key.width = unit(0.25,"cm"),legend.text =element_text(size=6),legend.key.height = unit(0.1,"cm"),legend.position = "bottom",legend.direction = "horizontal",legend.box="horizontal",legend.spacing.x = unit(0.08,"cm"),legend.box.margin=margin(-10, 0, 0, -100))+
  labs(col=expression("-log"[10]*italic("P")))+xlab("Counts")+
  scale_color_gradient(low="#F4B570", high = "#F06B37", limits=c(15,30))+guides(color=guide_colorbar(order=1),size=guide_legend(order=2))
ggsave("GO_down_24.pdf", height=60, units = "mm", width=80)

ggplot(data=GO_up_24[1:15,], aes(x=counts, y=reorder(Description,counts),size=ratio,color=LogP_pos))+geom_point()+
  scale_size(range = c(0.5, 3.5), name="Gene Ratio",breaks = c(0.2,0.3,0.4))+
  theme(axis.text.x=element_text(size=6.5,  family="sans",color="black"),axis.text.y=element_text(size=6.5,  family="sans",color="black"), 
        axis.title =element_text(size=6.5,  family="sans", color="black"),axis.line = element_line(colour = 'black', size = 0.46),axis.ticks.length=unit(.15, "cm"),
        axis.ticks =element_line(color="black",size=0.46),axis.title.y=element_blank(),panel.grid.major=element_line(size=0.46),legend.title = element_text(size=6.5),legend.key.width = unit(0.25,"cm"),legend.text =element_text(size=6),legend.key.height = unit(0.1,"cm"),legend.position = "bottom",legend.direction = "horizontal",legend.box="horizontal",legend.spacing.x = unit(0.08,"cm"),legend.box.margin=margin(-10, 0, 0, -100))+
  labs(col=expression("-log"[10]*italic("P")))+xlab("Counts")+
  scale_color_gradient(low="#F4B570", high = "#F06B37")+guides(color=guide_colorbar(order=1),size=guide_legend(order=2))
ggsave("GO_up_24.pdf", height=60, units = "mm", width=90)


