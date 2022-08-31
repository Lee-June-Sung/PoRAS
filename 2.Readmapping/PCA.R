#install.packages("devtools")
#install.packages("ggalt")
#install.packages('rlang')
#install.packages('digest')

library(devtools)
library(ggplot2)
library(RColorBrewer)


getwd()
setwd('D:/USER-DATA/Desktop/학사학위논문/Figure2/')
#setwd('D:/USER-DATA/Desktop/학사학위논문/Figure2')
#setwd("C:/Users/JangHG/Desktop")

x <- read.table('Figure2b_CDS (Total_FPKM).txt', row.names = 1, header = T)
group <- read.table("Figure2b_CDS_PCA_Group_name.txt",header=TRUE)


colors13 <- c(rgb(128,128,128, maxColorValue=255),    # 회
              rgb(237,0,0, maxColorValue=255),        # 빨
              rgb(0,84,255, maxColorValue=255),       # 파
              rgb(99,170,0, maxColorValue=255),       # 초
              rgb(255,187,0, maxColorValue=255),      # 노
              
              
              rgb(0,0,0, maxColorValue=255),          # 검
              rgb(255,112,18, maxColorValue=255),     # 주 
              rgb(113,18,255, maxColorValue=255),     # 남
              rgb(243,97,166, maxColorValue=255),     # 보
              rgb(43,165,186, maxColorValue=255))     #

log.x=log2(x+1)
xx = prcomp(t(log.x))
mPcv <- (xx$sdev^2/sum(xx$sdev^2))*100
names(mPcv) <- paste0("PC", 1:length(mPcv))
mDf <- data.frame(PC=1:10, var=mPcv[1:10])
#mScree <- ggplot(data=mDf, aes(x=PC, y=var)) +
#  geom_line(color="grey")+
#  geom_point() +
#  theme_bw(base_size = 20) +
#  scale_x_continuous(limits=c(1,10), breaks=1:10) +
#  ylim(c(100,150)) +
#  ggtitle("") + ylab("% explained variance")

mPca.df <- data.frame(Class=as.character(group$Class), xx$x)
mPca.df <- cbind(group$Group,mPca.df)
mPca.df <- cbind(group$Group2,mPca.df)

mPca.df

write.csv(mPca.df, file = 'test.csv')

(group$Group2)


levels(group$Group2)
levels(group$Group)
group$Group2 <-factor(group$Group2, levels=c("0h","6h","12h","1D","2D","3.5D","5D"))
group$Group <-factor(group$Group, levels=c("TDW","Pi"))


ggplot(mPca.df, aes(x=PC1,y=PC2, color=group$Group,shape=group$Group2))+ geom_point(size=3.5) +
  theme_bw(base_size = 12) + labs(color="Treatment", shape="Timepoint") +
  scale_shape_manual(values = c(3,1,2,0,16,17,15))+
  theme(axis.text=element_text(color = 'black', size=15),axis.title=element_text(size=15),legend.key.size = unit(0.5,"cm"),legend.key.height = unit(0.4,"cm"),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  coord_fixed(ratio = 1.4) + scale_color_manual(values = colors13) +
  stat_ellipse(aes(group = group$Group), size = 0.5, level = 0.95)+# + xlab('') + ylab('') 
  xlab(paste0("PC1 (", round(mPcv["PC1"], digits=2), " %)")) +#ylim(-120,120) +xlim(-100,160) + #cold만 xlim
  ylab(paste0("PC2 (", round(mPcv["PC2"], digits=2), " %)"))# +theme_bw(base_size = 15)
#




#res12 <- ggplot(mPca.df, aes(x=PC1,y=PC2,color=group$Group,shape=group$Group2)) + geom_point(size=3) +
#scale_shape_manual(values=1:length(unique(group$Class))) +
#  xlab(paste0("PC1 (", round(mPcv["PC1"], digits=2), " %)")) +
#  ylab(paste0("PC2 (", round(mPcv["PC2"], digits=2), " %)")) +
# theme_bw(base_size = 10) + ggtitle("") + labs(color="Treatment", shape="Timepoint") +
#theme(axis.text=element_text(size=15),axis.title=element_text(size=10),legend.key.size = unit(0.2,"cm"),legend.key.height = unit(0.2,"cm"),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
#scale_color_manual(values = colors13)  + coord_fixed(ratio = 1.4) #+  
#stat_ellipse(type = t)
#stat_ellipse(aes(x=PC1,y=PC2,color=factor(group$Group)))

#res12


#geom_encircle(aes(group = group$Group),linetype = 1, s_shape=1.2, expand=0.05, colour = group$Group)
#stat_ellipse(aes(x=PC1,y=PC2,color=factor(group$Group))) #,type = "norm", geom="polygon", level=0.8,alpha=0,size = 0.25)

#res12 + geom_encircle(data = group$Group, aes(x=PC1, y=PC2)) + 
# geom_encircle(data = group$Group, aes(x=PC1, y=PC2))

#res12 + geom_encircle(aes(group = group$Group))


#res12 + stat_ellipse(aes(x=PC1,y=PC2, color = factor(group$Group)))

#p <- ggord(res12, group$Group)

#ppi <- 300
#png("PC12A.png", width=5*ppi, height=5*ppi, res=ppi)
#plot(res12)
#dev.off()
#res12

#while (!is.null(dev.list()))  #dev.off()

