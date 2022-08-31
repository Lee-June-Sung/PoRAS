library(ggplot2)

getwd()                           # 본인의 경로 확인
setwd("D:/USER-DATA/Desktop/학사학위논문/")  # 경로 임의 설정

df <- read.table("Table/FIgure1c_CDS (Total_ratio).txt", header = T, row.names = 1, fill = T)

#barplot(prop.table(t(df), 2), las = 2, col = colorRampPalette(c("#FFA63B","#FFb08C","#FFE093","#263f44","#017991","#d9eeec"))(n=6))
barplot(prop.table(t(df[3:1]), 2), las = 2, col = colorRampPalette(c("#FFE093","#263f44","#017991","#d9eeec"))(n=4))
