use strict;

my ($fq1, $fq2);

open (DATA, "Read_statistics");
open (OUT, ">Read_statistics_split");
chomp(my $stLine = <DATA>);
print OUT $stLine."\t$stLine\t$stLine\n";
while ($stLine = <DATA>)
{
	chomp($stLine);
	if (($stLine =~ /1.fq/)||($stLine =~ /1.fastq/)||($stLine =~ /_1.cutadapt/)||($stLine =~ /1_non/))	
	{
		if ($stLine =~ /R1/)
		{
			$fq1 = $stLine;
		}
		elsif ($stLine =~ /R2/)
		{
			$fq2 = $stLine;
			my @list = split(/\t/, $fq1);
			my @list2 = split(/\t/, $fq2);
			my $fq3;
			my $calc = $list[1] + $list2[1];
			$fq3 = $calc;
			$calc = $list[2] + $list2[2]; 
			$fq3 = $fq3."\t$calc";
			$calc = ($list[3] + $list2[3])/2; 
			$fq3 = $fq3."\t$calc";
			$calc = ($list[4] + $list2[4])/2; 
			$fq3 = $fq3."\t$calc";
			$calc = ($list[5] + $list2[5])/2; 
			$fq3 = $fq3."\t$calc";
			$calc = $list[0];
			$calc =~ s/1.fq/Total/;
			$calc =~ s/1.fastq/Total/;
			$fq3 = $calc."\t$fq3";
			print OUT $fq1."\t$fq2\t$fq3\n";	
		}	
		$fq1 = $stLine;
	}
	elsif (($stLine =~ /2.fq/)||($stLine =~ /2.fastq/)||($stLine =~ /_2.cutadapt/)||($stLine =~ /2_non/))	
	{
		$fq2 = $stLine;
		my @list = split(/\t/, $fq1);
		my @list2 = split(/\t/, $fq2);
		my $fq3;
		my $calc = $list[1] + $list2[1];
		$fq3 = $calc;
		$calc = $list[2] + $list2[2]; 
		$fq3 = $fq3."\t$calc";
		$calc = ($list[3] + $list2[3])/2; 
		$fq3 = $fq3."\t$calc";
		$calc = ($list[4] + $list2[4])/2; 
		$fq3 = $fq3."\t$calc";
		$calc = ($list[5] + $list2[5])/2; 
		$fq3 = $fq3."\t$calc";
		$calc = $list[0];
		$calc =~ s/1.fq/Total/;
		$calc =~ s/1.fastq/Total/;
		$fq3 = $calc."\t$fq3";
		print OUT $fq1."\t$fq2\t$fq3\n";	
	}
}
