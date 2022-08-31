getwd()
setwd("D:/USER-DATA/Desktop/학사학위논문/Table/")

x <- read.table('FIgure1d_CDS (Total_FPKM).txt', row.names = 1, header = T)

boxplot(log2(x+1), las = 2, ylab = 'Log2 (FPKM+1)', cex=0.85, par(cex.axis=0.85, mar = c(14, 11, 7.5, 9)-7), ylim = c(0,16) ,col = c('#EAEAEA','#EAEAEA','#EAEAEA','#EAEAEA','#EAEAEA','#EAEAEA','#EAEAEA','#EAEAEA','#EAEAEA','#EAEAEA','#EAEAEA','#EAEAEA','#EAEAEA','#EAEAEA','#EAEAEA','#EAEAEA','#EAEAEA','#EAEAEA','#EAEAEA','#EAEAEA','#EAEAEA','#FF0000','#FF0000','#FF0000','#FF0000','#FF0000','#FF0000','#FF0000','#FF0000','#FF0000','#FF0000','#FF0000','#FF0000','#FF0000','#FF0000','#FF0000','#FF0000','#FF0000','#FF0000','#FF0000'))

# PDF로 11 x 7 inch 출력
