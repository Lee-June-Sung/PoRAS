setwd("D:\\Desktop\\DEG_Pi\\GO")
countdata <- read.table("Merged_BP.txt", header = T, sep = "\t")
library(ggplot2)
#y <-countdata[countdata$X.Log.FDR.>= 0,]
#countdata <- y
bubble <- ggplot(countdata, aes(x = reorder(Description,Sort), y = Group, color = X.LOG.P.value.)) + 
  geom_point(aes(size = numDEInCat)) + 
  scale_size(breaks = c(5,15,30), range = c(2,15)) +
  scale_color_gradientn(colors = c( "#90FFFF", "#050099")) +
  theme_bw() + 
  ylab("") +
  xlab("") + 
  labs(colour="-Log(FDR)", size="Count of Genes") +
  theme(strip.text = element_text(size = 15),legend.title = element_text(size = 15), axis.text = element_text(angle = 90, hjust = 1, vjust = 0.5,size=10),axis.title = element_text(size = 15),legend.text = element_text(size = 15) ,panel.border = element_rect(colour = "black"), panel.grid.major = element_line("#D8D8D8"), panel.grid.minor = element_line("white"))

bubble


