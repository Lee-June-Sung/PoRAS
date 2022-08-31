library(DESeq2)                                                               # DESeq2 library

getwd()                                                                       # current directory
setwd("D:/USER-DATA/Desktop/??????????????????/Figure3/")          # change directory

countdata <- read.table("CM-pi-6h_vs_CM-TDW-6h", header=TRUE, row.names=1) # Count file
#countdata <- read.table("CM-pi-12h_vs_CM-TDW-12h", header=TRUE, row.names=1) # Count file
#countdata <- read.table("CM-pi-1D_vs_CM-TDW-1D", header=TRUE, row.names=1) # Count file
#countdata <- read.table("CM-pi-2D_vs_CM-TDW-2D", header=TRUE, row.names=1) # Count file
#countdata <- read.table("CM-pi-3.5D_vs_CM-TDW-3.5D", header=TRUE, row.names=1) # Count file
#countdata <- read.table("CM-pi-5D_vs_CM-TDW-5D", header=TRUE, row.names=1) # Count file

countdata <- as.matrix(countdata)
countdata = countdata+1

(condition <- factor(c(rep("ctl", 3), rep("exp", 3))))
(coldata <- data.frame(row.names=colnames(countdata), condition))

dds <- DESeqDataSetFromMatrix(countData=round(countdata), colData=coldata, design=~condition)
dds <- DESeq(dds)
res <- results(dds)

table(res$log2FoldChange >= 2.0 & res$padj<0.05)
table(res$log2FoldChange <= -2.0 & res$padj<0.05)
res <- res[order(res$padj), ]

###############################

with(res, plot(log2FoldChange, -log10(padj), pch=19, cex=1.5, cex.lab=1.5, cex.axis=1.5, main="", xlim = c(-12,12), ylim = c(0,150)))

abline(h = -log10(0.05), col = "blue", lty = 2, lwd = 1)     # ???? ?????? ????
abline(v = c(-2,2), col = "blue", lty = 2, lwd = 1) # ???? ?????? ????

with(subset(res, -log10(padj) < -log10(0.05)), points(log2FoldChange, -log10(padj), pch=19, cex=1.5,col="darkgray")) # ???? ?????? ???? ?ƒÈ??? ???? ????
with(subset(res, log2FoldChange> -2 & log2FoldChange< 2 & -log10(padj)>-log10(0.05)), points(log2FoldChange, -log10(padj), pch=19,cex=1.5, col="darkgray")) # ???? ?????? ???? ?ƒÈ??? ???? ????
with(subset(res, log2FoldChange> 2 & -log10(padj)>-log10(0.05)), points(log2FoldChange, -log10(padj), pch=19, cex=1.5, col="#B42222"))        # up regulation = red
with(subset(res, log2FoldChange< -2 & -log10(padj)>-log10(0.05)), points(log2FoldChange, -log10(padj), pch=19, cex=1.5,col="#135AA1"))  # down regulation = blue 

