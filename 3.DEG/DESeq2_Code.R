
setwd("./DEG/")
library(DESeq2)

countdata <- read.table("CM-pi-12h_vs_CM-TDW-12h", header=TRUE, row.names=1)
countdata <- as.matrix(countdata)
head(countdata)
(condition <- factor(c(rep("ctl", 3), rep("exp", 3))))
(coldata <- data.frame(row.names=colnames(countdata), condition))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds <- DESeq(dds)
res <- results(dds)
table(res$padj<0.05)
res <- res[order(res$padj), ]
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head(resdata)
write.csv(resdata, file="DEG_result/CM-pi-12h_vs_CM-TDW-12h-results.csv")
resGA <- results(dds, lfcThreshold=0, altHypothesis="greaterAbs", alpha=0.05)
pdf("MA_plot/CM-pi-12h_vs_CM-TDW-12h.pdf",pointsize = 20)
drawLines <- function() abline(h=c(-1,1),col="dodgerblue",lwd=2)
plotMA(resGA, ylim=c(-5,5), main="CM-pi-12h_vs_CM-TDW-12h");drawLines()
dev.off()
a <- -log10(res$padj)
b <- a[is.finite(a)]
Up1 <- b[1]
Up2 <- b[1]*0.95
Up3 <- b[1]*0.9
Down1 <- b[1]
Down2 <- b[1]*0.95
Down3 <- b[1]*0.9
DownC1 <- length(row.names(subset(res, log2FoldChange< -1 & -log10(padj)>-log10(0.05))))
DownC2 <- length(row.names(subset(res, log2FoldChange< -2 & -log10(padj)>-log10(0.05))))
DownC3 <- length(row.names(subset(res, log2FoldChange< -4 & -log10(padj)>-log10(0.05))))
UpC1 <- length(row.names(subset(res, log2FoldChange> 1 & -log10(padj)>-log10(0.05))))
UpC2 <- length(row.names(subset(res, log2FoldChange> 2 & -log10(padj)>-log10(0.05))))
UpC3 <- length(row.names(subset(res, log2FoldChange> 4 & -log10(padj)>-log10(0.05))))
DownC1 = DownC1 - DownC2
DownC2 = DownC2 - DownC3
UpC1 = UpC1 - UpC2
UpC2 = UpC2 - UpC3
pdf("Volcano_plot/CM-pi-12h_vs_CM-TDW-12h.pdf",pointsize = 15)
with(res, plot(log2FoldChange, -log10(padj), pch=19, cex=1.5, cex.lab=1.5, cex.axis=1.5, main="", xlim = c(-10,10)))+
abline(h = -log10(0.05), col = "blue", lty = 2, lwd = 1)+     # 가로 경계선 출력
abline(v = c(-1,1), col = "blue", lty = 2, lwd = 1)+ # 세로 경계선 출력
abline(v = c(-2,2), col = "blue", lty = 2, lwd = 1)+ # 세로 경계선 출력
abline(v = c(-4,4), col = "blue", lty = 2, lwd = 1)+ # 세로 경계선 출력
with(subset(res, -log10(padj) < -log10(0.05)), points(log2FoldChange, -log10(padj), pch=19, cex=1.5,col="darkgray"))+
with(subset(res, log2FoldChange> -1 & log2FoldChange< 1 & -log10(padj)>-log10(0.05)), points(log2FoldChange, -log10(padj), pch=19,cex=1.5, col="darkgray")) +
with(subset(res, log2FoldChange> 1 & -log10(padj)>-log10(0.05)), points(log2FoldChange, -log10(padj), pch=19, cex=1.5, col="#FFC6C6"))+        # up regulation
with(subset(res, log2FoldChange> 2 & -log10(padj)>-log10(0.05)), points(log2FoldChange, -log10(padj), pch=19, cex=1.5, col="#FF7E7E"))+        # up regulation
with(subset(res, log2FoldChange> 4 & -log10(padj)>-log10(0.05)), points(log2FoldChange, -log10(padj), pch=19, cex=1.5, col="#FF0000"))+        # up regulation
with(subset(res, log2FoldChange< -1 & -log10(padj)>-log10(0.05)), points(log2FoldChange, -log10(padj), pch=19, cex=1.5,col="#ABF200"))+ # down regulation 
with(subset(res, log2FoldChange< -2 & -log10(padj)>-log10(0.05)), points(log2FoldChange, -log10(padj), pch=19, cex=1.5,col="#1DDB16"))+ # down regulation
with(subset(res, log2FoldChange< -4 & -log10(padj)>-log10(0.05)), points(log2FoldChange, -log10(padj), pch=19, cex=1.5,col="#008100"))+ # down regulation 
text(x=-9, y=Down1, label=DownC1, col="#ABF200")+ 
text(x=-9, y=Down2, label=DownC2, col="#1DDB16")+
text(x=-9, y=Down3, label=DownC3, col="#008100")+
text(x=9, y=Up1, label=UpC1, col="#FFC6C6")+
text(x=9, y=Up2, label=UpC2, col="#FF7E7E")+
text(x=9, y=Up3, label=UpC3, col="#FF0000")
dev.off()

countdata <- read.table("CM-pi-1D_vs_CM-TDW-1D", header=TRUE, row.names=1)
countdata <- as.matrix(countdata)
head(countdata)
(condition <- factor(c(rep("ctl", 3), rep("exp", 3))))
(coldata <- data.frame(row.names=colnames(countdata), condition))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds <- DESeq(dds)
res <- results(dds)
table(res$padj<0.05)
res <- res[order(res$padj), ]
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head(resdata)
write.csv(resdata, file="DEG_result/CM-pi-1D_vs_CM-TDW-1D-results.csv")
resGA <- results(dds, lfcThreshold=0, altHypothesis="greaterAbs", alpha=0.05)
pdf("MA_plot/CM-pi-1D_vs_CM-TDW-1D.pdf",pointsize = 20)
drawLines <- function() abline(h=c(-1,1),col="dodgerblue",lwd=2)
plotMA(resGA, ylim=c(-5,5), main="CM-pi-1D_vs_CM-TDW-1D");drawLines()
dev.off()
a <- -log10(res$padj)
b <- a[is.finite(a)]
Up1 <- b[1]
Up2 <- b[1]*0.95
Up3 <- b[1]*0.9
Down1 <- b[1]
Down2 <- b[1]*0.95
Down3 <- b[1]*0.9
DownC1 <- length(row.names(subset(res, log2FoldChange< -1 & -log10(padj)>-log10(0.05))))
DownC2 <- length(row.names(subset(res, log2FoldChange< -2 & -log10(padj)>-log10(0.05))))
DownC3 <- length(row.names(subset(res, log2FoldChange< -4 & -log10(padj)>-log10(0.05))))
UpC1 <- length(row.names(subset(res, log2FoldChange> 1 & -log10(padj)>-log10(0.05))))
UpC2 <- length(row.names(subset(res, log2FoldChange> 2 & -log10(padj)>-log10(0.05))))
UpC3 <- length(row.names(subset(res, log2FoldChange> 4 & -log10(padj)>-log10(0.05))))
DownC1 = DownC1 - DownC2
DownC2 = DownC2 - DownC3
UpC1 = UpC1 - UpC2
UpC2 = UpC2 - UpC3
pdf("Volcano_plot/CM-pi-1D_vs_CM-TDW-1D.pdf",pointsize = 15)
with(res, plot(log2FoldChange, -log10(padj), pch=19, cex=1.5, cex.lab=1.5, cex.axis=1.5, main="", xlim = c(-10,10)))+
abline(h = -log10(0.05), col = "blue", lty = 2, lwd = 1)+     # 가로 경계선 출력
abline(v = c(-1,1), col = "blue", lty = 2, lwd = 1)+ # 세로 경계선 출력
abline(v = c(-2,2), col = "blue", lty = 2, lwd = 1)+ # 세로 경계선 출력
abline(v = c(-4,4), col = "blue", lty = 2, lwd = 1)+ # 세로 경계선 출력
with(subset(res, -log10(padj) < -log10(0.05)), points(log2FoldChange, -log10(padj), pch=19, cex=1.5,col="darkgray"))+
with(subset(res, log2FoldChange> -1 & log2FoldChange< 1 & -log10(padj)>-log10(0.05)), points(log2FoldChange, -log10(padj), pch=19,cex=1.5, col="darkgray")) +
with(subset(res, log2FoldChange> 1 & -log10(padj)>-log10(0.05)), points(log2FoldChange, -log10(padj), pch=19, cex=1.5, col="#FFC6C6"))+        # up regulation
with(subset(res, log2FoldChange> 2 & -log10(padj)>-log10(0.05)), points(log2FoldChange, -log10(padj), pch=19, cex=1.5, col="#FF7E7E"))+        # up regulation
with(subset(res, log2FoldChange> 4 & -log10(padj)>-log10(0.05)), points(log2FoldChange, -log10(padj), pch=19, cex=1.5, col="#FF0000"))+        # up regulation
with(subset(res, log2FoldChange< -1 & -log10(padj)>-log10(0.05)), points(log2FoldChange, -log10(padj), pch=19, cex=1.5,col="#ABF200"))+ # down regulation 
with(subset(res, log2FoldChange< -2 & -log10(padj)>-log10(0.05)), points(log2FoldChange, -log10(padj), pch=19, cex=1.5,col="#1DDB16"))+ # down regulation
with(subset(res, log2FoldChange< -4 & -log10(padj)>-log10(0.05)), points(log2FoldChange, -log10(padj), pch=19, cex=1.5,col="#008100"))+ # down regulation 
text(x=-9, y=Down1, label=DownC1, col="#ABF200")+ 
text(x=-9, y=Down2, label=DownC2, col="#1DDB16")+
text(x=-9, y=Down3, label=DownC3, col="#008100")+
text(x=9, y=Up1, label=UpC1, col="#FFC6C6")+
text(x=9, y=Up2, label=UpC2, col="#FF7E7E")+
text(x=9, y=Up3, label=UpC3, col="#FF0000")
dev.off()

countdata <- read.table("CM-pi-2D_vs_CM-TDW-2D", header=TRUE, row.names=1)
countdata <- as.matrix(countdata)
head(countdata)
(condition <- factor(c(rep("ctl", 3), rep("exp", 3))))
(coldata <- data.frame(row.names=colnames(countdata), condition))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds <- DESeq(dds)
res <- results(dds)
table(res$padj<0.05)
res <- res[order(res$padj), ]
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head(resdata)
write.csv(resdata, file="DEG_result/CM-pi-2D_vs_CM-TDW-2D-results.csv")
resGA <- results(dds, lfcThreshold=0, altHypothesis="greaterAbs", alpha=0.05)
pdf("MA_plot/CM-pi-2D_vs_CM-TDW-2D.pdf",pointsize = 20)
drawLines <- function() abline(h=c(-1,1),col="dodgerblue",lwd=2)
plotMA(resGA, ylim=c(-5,5), main="CM-pi-2D_vs_CM-TDW-2D");drawLines()
dev.off()
a <- -log10(res$padj)
b <- a[is.finite(a)]
Up1 <- b[1]
Up2 <- b[1]*0.95
Up3 <- b[1]*0.9
Down1 <- b[1]
Down2 <- b[1]*0.95
Down3 <- b[1]*0.9
DownC1 <- length(row.names(subset(res, log2FoldChange< -1 & -log10(padj)>-log10(0.05))))
DownC2 <- length(row.names(subset(res, log2FoldChange< -2 & -log10(padj)>-log10(0.05))))
DownC3 <- length(row.names(subset(res, log2FoldChange< -4 & -log10(padj)>-log10(0.05))))
UpC1 <- length(row.names(subset(res, log2FoldChange> 1 & -log10(padj)>-log10(0.05))))
UpC2 <- length(row.names(subset(res, log2FoldChange> 2 & -log10(padj)>-log10(0.05))))
UpC3 <- length(row.names(subset(res, log2FoldChange> 4 & -log10(padj)>-log10(0.05))))
DownC1 = DownC1 - DownC2
DownC2 = DownC2 - DownC3
UpC1 = UpC1 - UpC2
UpC2 = UpC2 - UpC3
pdf("Volcano_plot/CM-pi-2D_vs_CM-TDW-2D.pdf",pointsize = 15)
with(res, plot(log2FoldChange, -log10(padj), pch=19, cex=1.5, cex.lab=1.5, cex.axis=1.5, main="", xlim = c(-10,10)))+
abline(h = -log10(0.05), col = "blue", lty = 2, lwd = 1)+     # 가로 경계선 출력
abline(v = c(-1,1), col = "blue", lty = 2, lwd = 1)+ # 세로 경계선 출력
abline(v = c(-2,2), col = "blue", lty = 2, lwd = 1)+ # 세로 경계선 출력
abline(v = c(-4,4), col = "blue", lty = 2, lwd = 1)+ # 세로 경계선 출력
with(subset(res, -log10(padj) < -log10(0.05)), points(log2FoldChange, -log10(padj), pch=19, cex=1.5,col="darkgray"))+
with(subset(res, log2FoldChange> -1 & log2FoldChange< 1 & -log10(padj)>-log10(0.05)), points(log2FoldChange, -log10(padj), pch=19,cex=1.5, col="darkgray")) +
with(subset(res, log2FoldChange> 1 & -log10(padj)>-log10(0.05)), points(log2FoldChange, -log10(padj), pch=19, cex=1.5, col="#FFC6C6"))+        # up regulation
with(subset(res, log2FoldChange> 2 & -log10(padj)>-log10(0.05)), points(log2FoldChange, -log10(padj), pch=19, cex=1.5, col="#FF7E7E"))+        # up regulation
with(subset(res, log2FoldChange> 4 & -log10(padj)>-log10(0.05)), points(log2FoldChange, -log10(padj), pch=19, cex=1.5, col="#FF0000"))+        # up regulation
with(subset(res, log2FoldChange< -1 & -log10(padj)>-log10(0.05)), points(log2FoldChange, -log10(padj), pch=19, cex=1.5,col="#ABF200"))+ # down regulation 
with(subset(res, log2FoldChange< -2 & -log10(padj)>-log10(0.05)), points(log2FoldChange, -log10(padj), pch=19, cex=1.5,col="#1DDB16"))+ # down regulation
with(subset(res, log2FoldChange< -4 & -log10(padj)>-log10(0.05)), points(log2FoldChange, -log10(padj), pch=19, cex=1.5,col="#008100"))+ # down regulation 
text(x=-9, y=Down1, label=DownC1, col="#ABF200")+ 
text(x=-9, y=Down2, label=DownC2, col="#1DDB16")+
text(x=-9, y=Down3, label=DownC3, col="#008100")+
text(x=9, y=Up1, label=UpC1, col="#FFC6C6")+
text(x=9, y=Up2, label=UpC2, col="#FF7E7E")+
text(x=9, y=Up3, label=UpC3, col="#FF0000")
dev.off()

countdata <- read.table("CM-pi-3.5D_vs_CM-TDW-3.5D", header=TRUE, row.names=1)
countdata <- as.matrix(countdata)
head(countdata)
(condition <- factor(c(rep("ctl", 3), rep("exp", 3))))
(coldata <- data.frame(row.names=colnames(countdata), condition))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds <- DESeq(dds)
res <- results(dds)
table(res$padj<0.05)
res <- res[order(res$padj), ]
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head(resdata)
write.csv(resdata, file="DEG_result/CM-pi-3.5D_vs_CM-TDW-3.5D-results.csv")
resGA <- results(dds, lfcThreshold=0, altHypothesis="greaterAbs", alpha=0.05)
pdf("MA_plot/CM-pi-3.5D_vs_CM-TDW-3.5D.pdf",pointsize = 20)
drawLines <- function() abline(h=c(-1,1),col="dodgerblue",lwd=2)
plotMA(resGA, ylim=c(-5,5), main="CM-pi-3.5D_vs_CM-TDW-3.5D");drawLines()
dev.off()
a <- -log10(res$padj)
b <- a[is.finite(a)]
Up1 <- b[1]
Up2 <- b[1]*0.95
Up3 <- b[1]*0.9
Down1 <- b[1]
Down2 <- b[1]*0.95
Down3 <- b[1]*0.9
DownC1 <- length(row.names(subset(res, log2FoldChange< -1 & -log10(padj)>-log10(0.05))))
DownC2 <- length(row.names(subset(res, log2FoldChange< -2 & -log10(padj)>-log10(0.05))))
DownC3 <- length(row.names(subset(res, log2FoldChange< -4 & -log10(padj)>-log10(0.05))))
UpC1 <- length(row.names(subset(res, log2FoldChange> 1 & -log10(padj)>-log10(0.05))))
UpC2 <- length(row.names(subset(res, log2FoldChange> 2 & -log10(padj)>-log10(0.05))))
UpC3 <- length(row.names(subset(res, log2FoldChange> 4 & -log10(padj)>-log10(0.05))))
DownC1 = DownC1 - DownC2
DownC2 = DownC2 - DownC3
UpC1 = UpC1 - UpC2
UpC2 = UpC2 - UpC3
pdf("Volcano_plot/CM-pi-3.5D_vs_CM-TDW-3.5D.pdf",pointsize = 15)
with(res, plot(log2FoldChange, -log10(padj), pch=19, cex=1.5, cex.lab=1.5, cex.axis=1.5, main="", xlim = c(-10,10)))+
abline(h = -log10(0.05), col = "blue", lty = 2, lwd = 1)+     # 가로 경계선 출력
abline(v = c(-1,1), col = "blue", lty = 2, lwd = 1)+ # 세로 경계선 출력
abline(v = c(-2,2), col = "blue", lty = 2, lwd = 1)+ # 세로 경계선 출력
abline(v = c(-4,4), col = "blue", lty = 2, lwd = 1)+ # 세로 경계선 출력
with(subset(res, -log10(padj) < -log10(0.05)), points(log2FoldChange, -log10(padj), pch=19, cex=1.5,col="darkgray"))+
with(subset(res, log2FoldChange> -1 & log2FoldChange< 1 & -log10(padj)>-log10(0.05)), points(log2FoldChange, -log10(padj), pch=19,cex=1.5, col="darkgray")) +
with(subset(res, log2FoldChange> 1 & -log10(padj)>-log10(0.05)), points(log2FoldChange, -log10(padj), pch=19, cex=1.5, col="#FFC6C6"))+        # up regulation
with(subset(res, log2FoldChange> 2 & -log10(padj)>-log10(0.05)), points(log2FoldChange, -log10(padj), pch=19, cex=1.5, col="#FF7E7E"))+        # up regulation
with(subset(res, log2FoldChange> 4 & -log10(padj)>-log10(0.05)), points(log2FoldChange, -log10(padj), pch=19, cex=1.5, col="#FF0000"))+        # up regulation
with(subset(res, log2FoldChange< -1 & -log10(padj)>-log10(0.05)), points(log2FoldChange, -log10(padj), pch=19, cex=1.5,col="#ABF200"))+ # down regulation 
with(subset(res, log2FoldChange< -2 & -log10(padj)>-log10(0.05)), points(log2FoldChange, -log10(padj), pch=19, cex=1.5,col="#1DDB16"))+ # down regulation
with(subset(res, log2FoldChange< -4 & -log10(padj)>-log10(0.05)), points(log2FoldChange, -log10(padj), pch=19, cex=1.5,col="#008100"))+ # down regulation 
text(x=-9, y=Down1, label=DownC1, col="#ABF200")+ 
text(x=-9, y=Down2, label=DownC2, col="#1DDB16")+
text(x=-9, y=Down3, label=DownC3, col="#008100")+
text(x=9, y=Up1, label=UpC1, col="#FFC6C6")+
text(x=9, y=Up2, label=UpC2, col="#FF7E7E")+
text(x=9, y=Up3, label=UpC3, col="#FF0000")
dev.off()

countdata <- read.table("CM-pi-5D_vs_CM-TDW-5D", header=TRUE, row.names=1)
countdata <- as.matrix(countdata)
head(countdata)
(condition <- factor(c(rep("ctl", 3), rep("exp", 3))))
(coldata <- data.frame(row.names=colnames(countdata), condition))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds <- DESeq(dds)
res <- results(dds)
table(res$padj<0.05)
res <- res[order(res$padj), ]
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head(resdata)
write.csv(resdata, file="DEG_result/CM-pi-5D_vs_CM-TDW-5D-results.csv")
resGA <- results(dds, lfcThreshold=0, altHypothesis="greaterAbs", alpha=0.05)
pdf("MA_plot/CM-pi-5D_vs_CM-TDW-5D.pdf",pointsize = 20)
drawLines <- function() abline(h=c(-1,1),col="dodgerblue",lwd=2)
plotMA(resGA, ylim=c(-5,5), main="CM-pi-5D_vs_CM-TDW-5D");drawLines()
dev.off()
a <- -log10(res$padj)
b <- a[is.finite(a)]
Up1 <- b[1]
Up2 <- b[1]*0.95
Up3 <- b[1]*0.9
Down1 <- b[1]
Down2 <- b[1]*0.95
Down3 <- b[1]*0.9
DownC1 <- length(row.names(subset(res, log2FoldChange< -1 & -log10(padj)>-log10(0.05))))
DownC2 <- length(row.names(subset(res, log2FoldChange< -2 & -log10(padj)>-log10(0.05))))
DownC3 <- length(row.names(subset(res, log2FoldChange< -4 & -log10(padj)>-log10(0.05))))
UpC1 <- length(row.names(subset(res, log2FoldChange> 1 & -log10(padj)>-log10(0.05))))
UpC2 <- length(row.names(subset(res, log2FoldChange> 2 & -log10(padj)>-log10(0.05))))
UpC3 <- length(row.names(subset(res, log2FoldChange> 4 & -log10(padj)>-log10(0.05))))
DownC1 = DownC1 - DownC2
DownC2 = DownC2 - DownC3
UpC1 = UpC1 - UpC2
UpC2 = UpC2 - UpC3
pdf("Volcano_plot/CM-pi-5D_vs_CM-TDW-5D.pdf",pointsize = 15)
with(res, plot(log2FoldChange, -log10(padj), pch=19, cex=1.5, cex.lab=1.5, cex.axis=1.5, main="", xlim = c(-10,10)))+
abline(h = -log10(0.05), col = "blue", lty = 2, lwd = 1)+     # 가로 경계선 출력
abline(v = c(-1,1), col = "blue", lty = 2, lwd = 1)+ # 세로 경계선 출력
abline(v = c(-2,2), col = "blue", lty = 2, lwd = 1)+ # 세로 경계선 출력
abline(v = c(-4,4), col = "blue", lty = 2, lwd = 1)+ # 세로 경계선 출력
with(subset(res, -log10(padj) < -log10(0.05)), points(log2FoldChange, -log10(padj), pch=19, cex=1.5,col="darkgray"))+
with(subset(res, log2FoldChange> -1 & log2FoldChange< 1 & -log10(padj)>-log10(0.05)), points(log2FoldChange, -log10(padj), pch=19,cex=1.5, col="darkgray")) +
with(subset(res, log2FoldChange> 1 & -log10(padj)>-log10(0.05)), points(log2FoldChange, -log10(padj), pch=19, cex=1.5, col="#FFC6C6"))+        # up regulation
with(subset(res, log2FoldChange> 2 & -log10(padj)>-log10(0.05)), points(log2FoldChange, -log10(padj), pch=19, cex=1.5, col="#FF7E7E"))+        # up regulation
with(subset(res, log2FoldChange> 4 & -log10(padj)>-log10(0.05)), points(log2FoldChange, -log10(padj), pch=19, cex=1.5, col="#FF0000"))+        # up regulation
with(subset(res, log2FoldChange< -1 & -log10(padj)>-log10(0.05)), points(log2FoldChange, -log10(padj), pch=19, cex=1.5,col="#ABF200"))+ # down regulation 
with(subset(res, log2FoldChange< -2 & -log10(padj)>-log10(0.05)), points(log2FoldChange, -log10(padj), pch=19, cex=1.5,col="#1DDB16"))+ # down regulation
with(subset(res, log2FoldChange< -4 & -log10(padj)>-log10(0.05)), points(log2FoldChange, -log10(padj), pch=19, cex=1.5,col="#008100"))+ # down regulation 
text(x=-9, y=Down1, label=DownC1, col="#ABF200")+ 
text(x=-9, y=Down2, label=DownC2, col="#1DDB16")+
text(x=-9, y=Down3, label=DownC3, col="#008100")+
text(x=9, y=Up1, label=UpC1, col="#FFC6C6")+
text(x=9, y=Up2, label=UpC2, col="#FF7E7E")+
text(x=9, y=Up3, label=UpC3, col="#FF0000")
dev.off()

countdata <- read.table("CM-pi-6h_vs_CM-TDW-6h", header=TRUE, row.names=1)
countdata <- as.matrix(countdata)
head(countdata)
(condition <- factor(c(rep("ctl", 3), rep("exp", 3))))
(coldata <- data.frame(row.names=colnames(countdata), condition))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds <- DESeq(dds)
res <- results(dds)
table(res$padj<0.05)
res <- res[order(res$padj), ]
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head(resdata)
write.csv(resdata, file="DEG_result/CM-pi-6h_vs_CM-TDW-6h-results.csv")
resGA <- results(dds, lfcThreshold=0, altHypothesis="greaterAbs", alpha=0.05)
pdf("MA_plot/CM-pi-6h_vs_CM-TDW-6h.pdf",pointsize = 20)
drawLines <- function() abline(h=c(-1,1),col="dodgerblue",lwd=2)
plotMA(resGA, ylim=c(-5,5), main="CM-pi-6h_vs_CM-TDW-6h");drawLines()
dev.off()
a <- -log10(res$padj)
b <- a[is.finite(a)]
Up1 <- b[1]
Up2 <- b[1]*0.95
Up3 <- b[1]*0.9
Down1 <- b[1]
Down2 <- b[1]*0.95
Down3 <- b[1]*0.9
DownC1 <- length(row.names(subset(res, log2FoldChange< -1 & -log10(padj)>-log10(0.05))))
DownC2 <- length(row.names(subset(res, log2FoldChange< -2 & -log10(padj)>-log10(0.05))))
DownC3 <- length(row.names(subset(res, log2FoldChange< -4 & -log10(padj)>-log10(0.05))))
UpC1 <- length(row.names(subset(res, log2FoldChange> 1 & -log10(padj)>-log10(0.05))))
UpC2 <- length(row.names(subset(res, log2FoldChange> 2 & -log10(padj)>-log10(0.05))))
UpC3 <- length(row.names(subset(res, log2FoldChange> 4 & -log10(padj)>-log10(0.05))))
DownC1 = DownC1 - DownC2
DownC2 = DownC2 - DownC3
UpC1 = UpC1 - UpC2
UpC2 = UpC2 - UpC3
pdf("Volcano_plot/CM-pi-6h_vs_CM-TDW-6h.pdf",pointsize = 15)
with(res, plot(log2FoldChange, -log10(padj), pch=19, cex=1.5, cex.lab=1.5, cex.axis=1.5, main="", xlim = c(-10,10)))+
abline(h = -log10(0.05), col = "blue", lty = 2, lwd = 1)+     # 가로 경계선 출력
abline(v = c(-1,1), col = "blue", lty = 2, lwd = 1)+ # 세로 경계선 출력
abline(v = c(-2,2), col = "blue", lty = 2, lwd = 1)+ # 세로 경계선 출력
abline(v = c(-4,4), col = "blue", lty = 2, lwd = 1)+ # 세로 경계선 출력
with(subset(res, -log10(padj) < -log10(0.05)), points(log2FoldChange, -log10(padj), pch=19, cex=1.5,col="darkgray"))+
with(subset(res, log2FoldChange> -1 & log2FoldChange< 1 & -log10(padj)>-log10(0.05)), points(log2FoldChange, -log10(padj), pch=19,cex=1.5, col="darkgray")) +
with(subset(res, log2FoldChange> 1 & -log10(padj)>-log10(0.05)), points(log2FoldChange, -log10(padj), pch=19, cex=1.5, col="#FFC6C6"))+        # up regulation
with(subset(res, log2FoldChange> 2 & -log10(padj)>-log10(0.05)), points(log2FoldChange, -log10(padj), pch=19, cex=1.5, col="#FF7E7E"))+        # up regulation
with(subset(res, log2FoldChange> 4 & -log10(padj)>-log10(0.05)), points(log2FoldChange, -log10(padj), pch=19, cex=1.5, col="#FF0000"))+        # up regulation
with(subset(res, log2FoldChange< -1 & -log10(padj)>-log10(0.05)), points(log2FoldChange, -log10(padj), pch=19, cex=1.5,col="#ABF200"))+ # down regulation 
with(subset(res, log2FoldChange< -2 & -log10(padj)>-log10(0.05)), points(log2FoldChange, -log10(padj), pch=19, cex=1.5,col="#1DDB16"))+ # down regulation
with(subset(res, log2FoldChange< -4 & -log10(padj)>-log10(0.05)), points(log2FoldChange, -log10(padj), pch=19, cex=1.5,col="#008100"))+ # down regulation 
text(x=-9, y=Down1, label=DownC1, col="#ABF200")+ 
text(x=-9, y=Down2, label=DownC2, col="#1DDB16")+
text(x=-9, y=Down3, label=DownC3, col="#008100")+
text(x=9, y=Up1, label=UpC1, col="#FFC6C6")+
text(x=9, y=Up2, label=UpC2, col="#FF7E7E")+
text(x=9, y=Up3, label=UpC3, col="#FF0000")
dev.off()
