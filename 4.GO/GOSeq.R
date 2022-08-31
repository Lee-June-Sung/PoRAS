
setwd("D:\\Desktop\\Reserch\\Dual_RNA-seq\\2022.02.09_DEG")
library(goseq)

GO = read.table("GO_list.txt", header=TRUE,sep = "\t")
KEGG = read.table("KEGG_list.txt", header=TRUE,sep = "\t")

Gene_ID = read.table("XCV_length.txt")
len = read.table("Lab-12h_filtered.txt")
Gene_ID$V3 <- 0
bias <- Gene_ID$V2
attr(bias, "names") <- Gene_ID$V1
DEGeneID <- Gene_ID$V3
attr(DEGeneID, "names") <- Gene_ID$V1
DEGeneID[len$V1] <- 1
bias <- as.numeric(bias)
pwf = nullp(DEGeneID, bias.data=bias)
pwf
GO.wall <- goseq(pwf, gene2cat=GO[,c(1,2)])
KEGG.wall <- goseq(pwf, gene2cat=KEGG[,c(1,2)])
write.csv(KEGG.wall,"Lab-12h_filtered_KEGG.csv")
write.csv(GO.wall,"Lab-12h_filtered_GO.csv")
dev.off()
DEGeneID = ""
len = ""
Gene_ID = ""
bias = ""
pwf = ""

GO = read.table("GO_list.txt", header=TRUE,sep = "\t")
KEGG = read.table("KEGG_list.txt", header=TRUE,sep = "\t")

Gene_ID = read.table("XCV_length.txt")
len = read.table("Lab-24h_filtered.txt")
Gene_ID$V3 <- 0
bias <- Gene_ID$V2
attr(bias, "names") <- Gene_ID$V1
DEGeneID <- Gene_ID$V3
attr(DEGeneID, "names") <- Gene_ID$V1
DEGeneID[len$V1] <- 1
bias <- as.numeric(bias)
pwf = nullp(DEGeneID, bias.data=bias)
pwf
GO.wall <- goseq(pwf, gene2cat=GO[,c(1,2)])
KEGG.wall <- goseq(pwf, gene2cat=KEGG[,c(1,2)])
write.csv(KEGG.wall,"Lab-24h_filtered_KEGG.csv")
write.csv(GO.wall,"Lab-24h_filtered_GO.csv")
dev.off()
DEGeneID = ""
len = ""
Gene_ID = ""
bias = ""
pwf = ""
GO = read.table("GO_list.txt", header=TRUE,sep = "\t")
KEGG = read.table("KEGG_list.txt", header=TRUE,sep = "\t")

Gene_ID = read.table("XCV_length.txt")
len = read.table("Lab-48h_filtered.txt")
Gene_ID$V3 <- 0
bias <- Gene_ID$V2
attr(bias, "names") <- Gene_ID$V1
DEGeneID <- Gene_ID$V3
attr(DEGeneID, "names") <- Gene_ID$V1
DEGeneID[len$V1] <- 1
bias <- as.numeric(bias)
pwf = nullp(DEGeneID, bias.data=bias)
pwf
GO.wall <- goseq(pwf, gene2cat=GO[,c(1,2)])
KEGG.wall <- goseq(pwf, gene2cat=KEGG[,c(1,2)])
write.csv(KEGG.wall,"Lab-48h_filtered_KEGG.csv")
write.csv(GO.wall,"Lab-48h_filtered_GO.csv")
dev.off()
DEGeneID = ""
len = ""
Gene_ID = ""
bias = ""
pwf = ""
GO = read.table("GO_list.txt", header=TRUE,sep = "\t")
KEGG = read.table("KEGG_list.txt", header=TRUE,sep = "\t")

Gene_ID = read.table("XCV_length.txt")
len = read.table("Tru-12h_filtered.txt")
Gene_ID$V3 <- 0
bias <- Gene_ID$V2
attr(bias, "names") <- Gene_ID$V1
DEGeneID <- Gene_ID$V3
attr(DEGeneID, "names") <- Gene_ID$V1
DEGeneID[len$V1] <- 1
bias <- as.numeric(bias)
pwf = nullp(DEGeneID, bias.data=bias)
pwf
GO.wall <- goseq(pwf, gene2cat=GO[,c(1,2)])
KEGG.wall <- goseq(pwf, gene2cat=KEGG[,c(1,2)])
write.csv(KEGG.wall,"Tru-12h_filtered_KEGG.csv")
write.csv(GO.wall,"Tru-12h_filtered_GO.csv")
dev.off()
DEGeneID = ""
len = ""
Gene_ID = ""
bias = ""
pwf = ""
GO = read.table("GO_list.txt", header=TRUE,sep = "\t")
KEGG = read.table("KEGG_list.txt", header=TRUE,sep = "\t")

Gene_ID = read.table("XCV_length.txt")
len = read.table("Tru-24h_filtered.txt")
Gene_ID$V3 <- 0
bias <- Gene_ID$V2
attr(bias, "names") <- Gene_ID$V1
DEGeneID <- Gene_ID$V3
attr(DEGeneID, "names") <- Gene_ID$V1
DEGeneID[len$V1] <- 1
bias <- as.numeric(bias)
pwf = nullp(DEGeneID, bias.data=bias)
pwf
GO.wall <- goseq(pwf, gene2cat=GO[,c(1,2)])
KEGG.wall <- goseq(pwf, gene2cat=KEGG[,c(1,2)])
write.csv(KEGG.wall,"Tru-24h_filtered_KEGG.csv")
write.csv(GO.wall,"Tru-24h_filtered_GO.csv")
dev.off()
DEGeneID = ""
len = ""
Gene_ID = ""
bias = ""
pwf = ""
GO = read.table("GO_list.txt", header=TRUE,sep = "\t")
KEGG = read.table("KEGG_list.txt", header=TRUE,sep = "\t")

Gene_ID = read.table("XCV_length.txt")
len = read.table("Tru-48h_filtered.txt")
Gene_ID$V3 <- 0
bias <- Gene_ID$V2
attr(bias, "names") <- Gene_ID$V1
DEGeneID <- Gene_ID$V3
attr(DEGeneID, "names") <- Gene_ID$V1
DEGeneID[len$V1] <- 1
bias <- as.numeric(bias)
pwf = nullp(DEGeneID, bias.data=bias)
pwf
GO.wall <- goseq(pwf, gene2cat=GO[,c(1,2)])
KEGG.wall <- goseq(pwf, gene2cat=KEGG[,c(1,2)])
write.csv(KEGG.wall,"Tru-48h_filtered_KEGG.csv")
write.csv(GO.wall,"Tru-48h_filtered_GO.csv")
dev.off()
DEGeneID = ""
len = ""
Gene_ID = ""
bias = ""
pwf = ""

Gene_ID = read.table("XCV_length.txt")
len = read.table("Lab.txt")
Gene_ID$V3 <- 0
bias <- Gene_ID$V2
attr(bias, "names") <- Gene_ID$V1
DEGeneID <- Gene_ID$V3
attr(DEGeneID, "names") <- Gene_ID$V1
DEGeneID[len$V1] <- 1
bias <- as.numeric(bias)
pwf = nullp(DEGeneID, bias.data=bias)
pwf
GO.wall <- goseq(pwf, gene2cat=GO[,c(1,2)])
KEGG.wall <- goseq(pwf, gene2cat=KEGG[,c(1,2)])
write.csv(KEGG.wall,"Lab_KEGG.csv")
write.csv(GO.wall,"Lab_GO.csv")
dev.off()
DEGeneID = ""
len = ""
Gene_ID = ""
bias = ""
pwf = ""

Gene_ID = read.table("XCV_length.txt")
len = read.table("Tru.txt")
Gene_ID$V3 <- 0
bias <- Gene_ID$V2
attr(bias, "names") <- Gene_ID$V1
DEGeneID <- Gene_ID$V3
attr(DEGeneID, "names") <- Gene_ID$V1
DEGeneID[len$V1] <- 1
bias <- as.numeric(bias)
pwf = nullp(DEGeneID, bias.data=bias)
pwf
GO.wall <- goseq(pwf, gene2cat=GO[,c(1,2)])
KEGG.wall <- goseq(pwf, gene2cat=KEGG[,c(1,2)])
write.csv(KEGG.wall,"Tru_KEGG.csv")
write.csv(GO.wall,"Tru_GO.csv")
dev.off()
DEGeneID = ""
len = ""
Gene_ID = ""
bias = ""
pwf = ""

