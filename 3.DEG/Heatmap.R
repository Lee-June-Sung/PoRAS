# ComplexHeatmap install
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

#BiocManager::install("ComplexHeatmap")
########################
#install.packages('magick')

#library('magick')
library(ComplexHeatmap)
library('dplyr')
library(circlize)

getwd()
setwd("D:/USER-DATA/Desktop/")          # change directory
#setwd("D:/01.업무용/02. 과제 관련/00_Dual RNA-seq 프로젝트/Dual RNA-seq 미팅/21.08.17/")

#data1 = read.table("common_fpkm.txt", header = TRUE, sep = '\t', row.names="Gene_ID")    # common
#data1 = read.table("enriched_fpkm.txt", header = TRUE, sep = '\t', row.names="Gene_ID")     # enriched
data1 = read.table("log2FC_2_FPKM.txt", header = TRUE, sep = '\t', row.names="Gene_ID")     # conventional


#data1 = read.table("Heatmap/Total_Lab_vs_Tru/common_gene_FPKM.txt", header = TRUE, sep = '\t', row.names="Gene_ID")
#data1 = read.table("Heatmap/Total_Lab_vs_Tru/con_specific_gene_FPKM.txt", header = TRUE, sep = '\t', row.names="Gene_ID")


#data1 = read.table("Total_Lab24_pval-sorted_SNP_GO.txt", header = TRUE, sep = '\t', row.names = 'term')
#data1 = read.table("Total_Lab48_pval-sorted_SNP_GO.txt", header = TRUE, sep = '\t', row.names = 'term')

#data1 = read.table("Total_Tru12_pval-sorted_SNP_GO.txt", header = TRUE, sep = '\t', row.names = 'term')
#data1 = read.table("Total_Tru24_pval-sorted_SNP_GO.txt", header = TRUE, sep = '\t', row.names = 'term')
#data1 = read.table("Total_Tru48_pval-sorted_SNP_GO.txt", header = TRUE, sep = '\t', row.names = 'term')


#data1 = read.table("Total_Lab12_pval-sorted_SNP_GO.txt", header=TRUE, sep="\t", row.names="terms")
x <- data.matrix(data1)
x <- log2(x+1)
#x <- -log10(x)
#x <- log2(x+1)


#data2 = read.table("Heatmap/Total_Lab_vs_Tru/con_specific_gene_FC.txt", header=TRUE, sep="\t", row.names="Gene_ID")
#y <- data.matrix(data2)

data2 = read.table("log2FC_2_FC.txt", header=TRUE, sep="\t", row.names="Gene_ID")
y <- data.matrix(data2)

data3 = read.table("log2FC_2_FDR.txt", header=TRUE, sep="\t", row.names="Gene_ID")   # conventional
z <- -log10(data.matrix(data3))

#data3 = read.table("common_FDR.txt", header=TRUE, sep="\t", row.names="Gene_ID")   # common
#data3 = read.table("enriched_FDR.txt", header=TRUE, sep="\t", row.names="Gene_ID")    # enriched

# FDR
#data3 = read.table("conventional_FDR.txt", header=TRUE, sep="\t", row.names="Gene_ID")   # conventional
#z <- -log10(data.matrix(data3))

#data4 = read.table("common_fpkm_ratio_plus1.txt", header=TRUE, sep="\t", row.names="Gene_ID")    # common
#data4 = read.table("enriched_fpkm_ratio_plus1.txt", header=TRUE, sep="\t", row.names="Gene_ID")    # enriched

# FC
#data4 = read.table("conventional_fpkm_ratio_plus1.txt", header=TRUE, sep="\t", row.names="Gene_ID")    # conventional
#a <- -log2(data.matrix(data4))

#-log10(0.05)
#-log10(0.005)
#-log10(0.0005)
#-log10(0.00005)
#-log10(0.000005)

#Heatmap(subset(x), name = "-log10 pvalue", km = 1, col = colorRamp2(c(1.30102,1.30103,2.30103,3.30103,4.30103,5.30103), c("black","darkgreen","Yellow","Orange","red","#CE1212")), cluster_columns = F, cluster_rows =  F,
#        show_row_names = T, show_column_names = F) +


Heatmap(x, name = "Log2 FPKM", km = 1, heatmap_legend_param = list(at=c(0,1,2,3,4,5), color_bar = "continuous"), 
        col = colorRamp2(c(0,1,2,3,4,5), c("black","darkgreen","Yellow","Orange","red","darkred")),   # heatmap1
        cluster_columns = F, cluster_rows =  T, show_row_names = F, show_column_names = F) +
Heatmap(y, name = "Log2 FC", heatmap_legend_param = list(at=c(-4,-2,0,2,4), color_bar = "continuous"), col = circlize::colorRamp2(c(-4,-2,0,2,4),c("#4B96AA","#3C5064","black","#969628","#FAE650")), show_row_names = F, show_column_names = F, cluster_columns = F,cluster_rows =  F) +
Heatmap(z, name = "- Log10 FDR", heatmap_legend_param = list(at=c(0,3,6,9), color_bar = "continuous"), col = circlize::colorRamp2(c(0,3,6,9),c("black","#5c8d89","#74b49b","#a7d7c5")), show_row_names = F, show_column_names = F, cluster_columns = F,cluster_rows =  F) 

#Heatmap(x, name = "Log2 FPKM", km = 1, heatmap_legend_param = list(at=c(0,1.5,3,4.5,6,7.5), color_bar = "continuous"), col = colorRamp2(c(0,1.5,3,4.5,6,7.5), c("black","darkgreen","Yellow","Orange","red","darkred")),   # heatmap1
        #Heatmap(x, name = "Log2 FPKM", km = 1, heatmap_legend_param = list(at=c(0,2,4,6,8,10), color_bar = "continuous"), col = colorRamp2(c(0,2,4,6,8,10), c("black","darkgreen","Yellow","Orange","red","darkred")),   # heatmap2
        #Heatmap(x, name = "Log2 FPKM", km = 1, heatmap_legend_param = list(at=c(0,2.5,5,7.5,10,12.5), color_bar = "continuous"), col = colorRamp2(c(0,2.5,5,7.5,10,12.5), c("black","darkgreen","Yellow","Orange","red","darkred")),   # heatmap3
        #Heatmap(x, name = "Log2 FPKM", km = 1, col = colorRamp2(c(0,2,4,6,8,10), c("black","darkgreen","Yellow","Orange","red","darkred")), 
#        cluster_columns = F, cluster_rows =  T, show_row_names = F, show_column_names = F)
  #Heatmap(subset(x), name = "Log2 FP`KM", km = 1, col = colorRamp2(c(0,2,4,6,8,10), c("black","darkgreen","Yellow","Orange","red","darkred")), cluster_columns = F, cluster_rows =  F,
  #Heatmap(subset(x), name = "Log2 FPKM", km = 1, col = colorRamp2(c(0,2,4,6,8), c("black","darkgreen","Yellow","Orange","red")), cluster_columns = F, cluster_rows =  T,
  
  
  #Heatmap(y, name = "Log2 FC", heatmap_legend_param = list(at=c(-4,-2,0,2,4), color_bar = "continuous"), col = circlize::colorRamp2(c(-4,-2,0,2,4),c("#4B96AA","#3C5064","black","#969628","#FAE650")), show_row_names = F, show_column_names = F, cluster_columns = F,cluster_rows =  F) +
  #Heatmap(a, name = "Log2 FPKM ratio", heatmap_legend_param = list(at=c(0,1,2), color_bar = "continuous"), col = circlize::colorRamp2(c(0,1,2),c("black","#969628","#FAE650")), show_row_names = F, show_column_names = F, cluster_columns = F,cluster_rows =  F) +

# 이거 사용
#  ##Heatmap(a, name = "Log2 FPKM ratio", heatmap_legend_param = list(at=c(-1,0,1), color_bar = "continuous"), col = circlize::colorRamp2(c(-1,0,1),c("#4678B7","#FEFEBE","#DA372B")), show_row_names = F, show_column_names = F, cluster_columns = F,cluster_rows =  F) +
#  ##Heatmap(z, name = "- Log10 FDR", heatmap_legend_param = list(at=c(0,3,6,9), color_bar = "continuous"), col = circlize::colorRamp2(c(0,3,6,9),c("black","#5c8d89","#74b49b","#a7d7c5")), show_row_names = F, show_column_names = F, cluster_columns = F,cluster_rows =  F) 

# color_bar = "continuous" 대신 “discrete” 사용하면 그라데이션이 아니라 구간으로 나뉘게 됨



