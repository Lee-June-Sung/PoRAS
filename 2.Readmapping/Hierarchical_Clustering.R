source("http://www.bioconductor.org/biocLite.R")
biocLite("tweeDEseq")
install.packages("gplots")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("tweeDEseq")


getwd()
#setwd("D:/USER-DATA/Desktop/Signal molecule")
#setwd("E:/Dropbox (분자방)/분자방/2-1.전략_RNAseq/2. Dual RNA-seq_opt/Dual RNA-seq Figure/Hierarchical Clustering")
setwd("D:/USER-DATA/Desktop/학사학위논문/Figure2/")


library(DESeq2)
library(edgeR)
library(tweeDEseq)
library("gplots")

count.data <- read.table("FIgure2a_CDS (Total_FPKM).txt",sep="\t",header=T,row.names = 1)
sample.information <- read.table("FIgure2a_CDS (Total_FPKM).txt",sep="\t",header=T,row.names = 1)
count.data <- count.data+1
numberRead = 15
numberSample = 5
keep <- rowSums(count.data >= numberRead) > numberSample
filter.data <- count.data[keep,]
filter.data  <- normalizeCounts(filter.data ,group=sample.information$Group, method="TMM")
input.data <- log2(filter.data+1)
hdata = hclust(dist(t(input.data)),method = "com")
plot(hdata, hang = -1)
hdata = hclust(dist(t(input.data)),method = "ave")
plot(hdata, hang = -1)
hdata = hclust(dist(t(input.data)),method = "ward.D2")
plot(hdata, hang = -1)
dd=as.dist(1-cor(input.data))
plot(hclust(dd, method="com"))
cor.data = cor(input.data)
yb <- colorRampPalette(c("#FFFFFF","#FFFF24","#FF0000"))
heatmap(as.matrix(cor.data),margins = c(8, 8))
heatmap.2(as.matrix(cor.data),col=yb, trace = "none", margins = c(8,10),density.info = "none",key = F)
?heatmap.2

