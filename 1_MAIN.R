library("genefilter")
library("ggrepel")
library('ggupset')
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("plyr")
library("dplyr")
library("reshape")
library("biomaRt")
library("grid")
library("gridExtra")
library("readxl")
library("stringr")
library("reshape2")
library('clusterProfiler')
library('enrichplot')
library('org.Mm.eg.db')
library('ggnewscale')
library('devtools')
library('ComplexHeatmap')
library('tidyverse')
library("vsn")
library("R.utils")
library("optparse")
library('circlize')



library("ggplot2")
theme = theme(panel.background = element_blank(),
              panel.border=element_rect(fill=NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              strip.background=element_blank(),
              axis.text.x=element_text(colour="black",size=10),
              axis.text.y=element_text(colour="black",size=10),
              axis.ticks=element_line(colour="black"),
              plot.margin=unit(c(1,1,1,1),"line"))


#######################  Main script for all RNA-seq related analysis 

cts = read.csv('data/RSEMCountTable.tab',sep = '\t',header = TRUE)
cts = unique(cts)
row.names(cts) = cts$Gene
cts$Gene = NULL
cts$X = NULL
sample = colnames(cts)
coldata = as.data.frame(sample)
coldata$experiment = sapply(strsplit(as.character(coldata$sample), "_"), "[[", 1 )
coldata$tissue = sapply(strsplit(as.character(coldata$sample), "_"), "[[", 2 )
coldata$diet = sapply(strsplit(as.character(coldata$sample), "_"), "[[", 3 )
coldata$group = paste(coldata$experiment, coldata$tissue, coldata$diet,sep = '_')
coldata$tissue = gsub("Adj", "", coldata$tissue)
row.names(coldata) = coldata$sample

head(cts,n=2)
head(coldata,n=2)

#######################  Load DESeq2 object

dds = DESeqDataSetFromMatrix(countData = cts,colData = coldata,design= ~ group)

#######################  PCA plot - Figure 2C (write to output/Figures)

rld = rlog(dds, blind=FALSE)
pcaData = plotPCA(rld, intgroup=c('group','diet'),returnData=TRUE,ntop=99999999999999)
#pcaData$group = factor(pcaData$group,levels=level_order,ordered = TRUE)
unique(pcaData$group)
percentVar = round(100 * attr(pcaData, "percentVar"))
g = ggplot(pcaData, aes(PC1, PC2))
g = g + geom_point(aes(shape=group,color=diet),size=6)
g = g + scale_shape_manual(values = c(19,19,15,1,1))
g = g + scale_color_manual(values = c('blue','red'))
g = g + ggtitle("Principal component analysis")
g = g + xlab(paste0("PC1 (",percentVar[1],"% variance)"))
g = g + ylab(paste0("PC2 (",percentVar[2],"% variance)"))
g = g + coord_fixed() 
g = g + theme 
g = g + ggtitle('')
g = g + theme(text = element_text(size=15),
              axis.text.x = element_text(size = 20),
              axis.text.y = element_text(size = 20),  
              axis.title.x = element_text(size = 20),
              axis.title.y = element_text(size = 20))
g
pdf('Figure1/Figure2C.pdf')
print(g)
dev.off()

#######################  Perform differential gene expression

comps = c('nonDEN_Liver_WD_vs_nonDEN_Liver_CD',
          'DEN_Liver_CD_vs_nonDEN_Liver_CD',
          'DEN_AdjLiver_WD_vs_nonDEN_Liver_CD',
          'DEN_Tumour_WD_vs_nonDEN_Liver_CD',
          'DEN_AdjLiver_WD_vs_nonDEN_Liver_WD',
          'DEN_Liver_CD_vs_nonDEN_Liver_WD',
          'DEN_Tumour_WD_vs_nonDEN_Liver_WD',
          'DEN_AdjLiver_WD_vs_DEN_Liver_CD',
          'DEN_Tumour_WD_vs_DEN_AdjLiver_WD',
          'DEN_Tumour_WD_vs_DEN_Liver_CD')      # manually specified to ensure correct baseling for each DGE

counter = 0
for (comp in comps){
  counter = counter + 1
  
  #perform DESeq2
  ref = strsplit(comp,'_vs_')[[1]][2]
  dds = DESeqDataSetFromMatrix(countData = cts,colData = coldata,design= ~ group)
  dds$group = relevel(dds$group,ref = ref)
  dds = DESeq(dds)
  res = lfcShrink(dds,coef=paste('group',comp,sep='_'),type='apeglm')
  
  if(counter == 1) {
    genes = row.names(res)
    resTable = as.data.frame(genes)
  }
  
  #add in all the columns
  player = paste('contrast',counter,'lg2BaseMean',comp,sep='_')
  resTable[player] = res$baseMean
  
  player = paste('contrast',counter,'lg10p',comp,sep='_')
  resTable[player] = res$pvalue
  
  player = paste('contrast',counter,'padj',comp,sep='_')
  resTable[player] = res$padj
  
  player = paste('contrast',counter,'logFC',comp,sep='_')
  resTable[player] = res$log2FoldChange
}
# table with DEGs for each contrast
write.csv(resTable,'data/esLFC.csv')


#######################  get enriched GO Biological Processes  - Figure 2E (write to output/Figures)

combinedDotPlot = function(resLFC,homology,contrasts){
  gseDat = list()
  for (contrast in contrasts){
    search = paste('contrast_',contrast,'_',sep='')
    resLFC$entrez = homologyMap$NCBI.gene..formerly.Entrezgene..ID[match(resLFC$genes,homologyMap$Gene.stable.ID)]
    df = resLFC[ , grepl( search , names( resLFC ) ) ]
    parts = str_split(colnames(df)[1],'lg2BaseMean_')
    title = parts[[1]][2]
    print(title)
    df = data.frame(entrez = resLFC$entrez,lfc = df[,4],padj=df[,3])
    df = df[!is.na(df$entrez),]
    df = subset(df,abs(df$lfc) > 1 & df$padj <= 0.05)
    de = unique(df$entrez)
    gseDat[[title]] = de
  }
  ck = compareCluster(geneCluster = gseDat, 
                      fun = "enrichGO",
                      pAdjustMethod = "BH",
                      OrgDb = "org.Mm.eg.db", 
                      ont = "bp", 
                      readable = TRUE, 
                      qvalueCutoff  = 0.05)
  return(ck)
}

resLFC = read.csv('data/resLFC.csv',header = TRUE,sep = ',')
resLFC$X = NULL
homologyMap = read.csv('data/homologyMap.tab',header = TRUE, sep = '\t')
contrasts = c(1,2,3,4,5,6,7,8,9,10)
p = combinedDotPlot(resLFC,homology,contrasts)

level_order = c('DEN_Tumour_WD','DEN_AdjLiver_WD','DEN_Liver_CD','nonDEN_Liver_WD','nonDEN_Liver_CD')
comps = c('DEN_Tumour_WD_vs_DEN_AdjLiver_WD','DEN_Tumour_WD_vs_DEN_Liver_CD','DEN_Tumour_WD_vs_nonDEN_Liver_WD',
          'DEN_Tumour_WD_vs_nonDEN_Liver_CD','DEN_AdjLiver_WD_vs_DEN_Liver_CD','DEN_AdjLiver_WD_vs_nonDEN_Liver_WD',
          'DEN_AdjLiver_WD_vs_nonDEN_Liver_CD','DEN_Liver_CD_vs_nonDEN_Liver_WD','DEN_Liver_CD_vs_nonDEN_Liver_CD',
          'nonDEN_Liver_WD_vs_nonDEN_Liver_CD')

comps = rev(comps)
q = p
q = dotplot(q)
dat = q$data
names(dat)[names(dat) == "p.adjust"] <- "padj"
dat = na.omit(dat)
dat$Cluster = gsub("\n.*","",dat$Cluster)
dat$geneID = NULL
dat$Cluster = factor(dat$Cluster,levels=comps,ordered = TRUE)
dat = subset(dat,dat$padj <= 0.001) 
processOrder = read.csv('data/processOrder.csv',sep=',',header=T)
dat$Description = factor(dat$Description,levels=rev(processOrder$Process),ordered = T)

g = ggplot(data=dat,aes(x=Cluster,y=Description))
g = g + geom_point(aes(color=padj,size=GeneRatio))
g = g + theme
#g = g + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
g = g + theme(text = element_text(size=12),
              axis.text.x = element_text(size = 12),
              axis.text.y = element_text(size = 12),  
              axis.title.x = element_text(size = 12),
              axis.title.y = element_text(size = 12))
g = g + axis_combmatrix(sep = "_vs_",levels = level_order)
g = g + scale_color_gradient(low="red", high="darkblue")
g = g + ylab('')
g = g + xlab('DGE Experiment')
g
pdf('Figure1/Figure2D.pdf')
print(g)
dev.off()








