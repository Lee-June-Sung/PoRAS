library(ggplot2)

getwd()                           # ������ ��� Ȯ��
setwd("D:/USER-DATA/Desktop/�л���������/")  # ��� ���� ����

df <- read.table("Table/FIgure1c_CDS (Total_ratio).txt", header = T, row.names = 1, fill = T)

#barplot(prop.table(t(df), 2), las = 2, col = colorRampPalette(c("#FFA63B","#FFb08C","#FFE093","#263f44","#017991","#d9eeec"))(n=6))
barplot(prop.table(t(df[3:1]), 2), las = 2, col = colorRampPalette(c("#FFE093","#263f44","#017991","#d9eeec"))(n=4))