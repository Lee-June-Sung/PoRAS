use strict;

my $Set = $ARGV[0];

my (%bank, @DEG);

open (DATA, $Set);
while (my $stLine = <DATA>)
{
	chomp($stLine);
	$stLine =~ s///;
	my @list = split(/\t/, $stLine);
	my $set = $list[1]."_vs_".$list[0];
	push @DEG, $set;
}
close DATA;
system("mkdir DEG");
system("mkdir DEG/DEG_result");
system("mkdir DEG/MA_plot");
system("mkdir DEG/Volcano_plot");

open (OUT2, ">DESeq2_Code.txt");

print OUT2 '
setwd("./DEG/")
library(DESeq2)
';
for (my $i = 0; $i < @DEG; $i++)
{
	my @list = split(/_vs_/, $DEG[$i]);
	my $out = "DEG/$list[0]\_vs\_$list[1]";
	open (OUT, ">$out");
	open (DATA, "Paired_Read_count.txt");
	my $C = 0;
	my $T = 0;
	my @line1 = ();
	my @line2 = ();
	my $check4 = "OK";
	my $check5 = "OK";
	my @name1 = ();
	my @name2 = ();
	while (my $stLine = <DATA>)
	{
		chomp($stLine);
		$stLine =~ s///;
		my @list2 = split(/\t/, $stLine);
		my $check = "NO";
		my $check2 = "";
		my $check3 = (scalar @list2 -1)/3;
		for (my $j = 1; $j < @list2; $j++)
		{
			if ($list2[$j] >= 0)
			{
				$check2++;
			}
		}
#		if ($check2 > $check3)
#		{
			$check = "OK";
#		}
		if ($stLine =~ /Gene_ID/)
		{
			for (my $j = 0; $j < @list2; $j++)
			{
				if ($list2[$j] =~ /$list[0]/)
				{
					push @name1, $list2[$j];
					push @line1, $j;
				}
				if ($list2[$j] =~ /$list[1]/)
				{
					push @name2, $list2[$j];
					push @line2, $j;
				}
			}
		}
		if ($check eq "OK")
		{
			if ($check4 eq "OK")
			{
#				print OUT "\t";
#				print OUT $_."\t" foreach (@name2);
				$check4 = "NO";
			}
			if ($check5 eq "OK")
			{
#				print OUT $_."\t" foreach (@name1);
#				print OUT "\n";
				$check5 = "NO";
			}

			print OUT "$list2[0]";
			for (my $j = 0; $j < @line2; $j++)
			{
				print OUT "\t$list2[$line2[$j]]";
			}	
			for (my $j = 0; $j < @line1; $j++)
			{
				print OUT "\t$list2[$line1[$j]]";
			}	
			print OUT "\n";
		}
	}
	$C = scalar @line2;
	$T = scalar @line1;
	close DATA;
	print OUT2 '
countdata <- read.table("'."$list[0]_vs_$list[1]".'", header=TRUE, row.names=1)
countdata <- as.matrix(countdata)
head(countdata)
(condition <- factor(c(rep("ctl", '."$C".'), rep("exp", '."$T".'))))
(coldata <- data.frame(row.names=colnames(countdata), condition))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds <- DESeq(dds)
res <- results(dds)
table(res$padj<0.05)
res <- res[order(res$padj), ]
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head(resdata)
write.csv(resdata, file="'."DEG_result/$list[0]_vs_$list[1]".'-results.csv")
resGA <- results(dds, lfcThreshold=0, altHypothesis="greaterAbs", alpha=0.05)
pdf("'."MA_plot/$list[0]_vs_$list[1].pdf".'",pointsize = 20)
drawLines <- function() abline(h=c(-1,1),col="dodgerblue",lwd=2)
plotMA(resGA, ylim=c(-5,5), main="'."$list[0]_vs_$list[1]".'");drawLines()
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
pdf("'."Volcano_plot/$list[0]_vs_$list[1].pdf".'",pointsize = 15)
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
';


}

