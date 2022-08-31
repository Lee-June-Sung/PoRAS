use strict;

my $Model = $ARGV[0]; # input : Gene Model or Genome
my @fastq = 0;

my ($fq1, $fq2);
open (OUT, ">1.Hisat2.txt");

@fastq = glob("*.fq");
if (scalar @fastq == 0)
{
	@fastq = glob("*.fq.gz");
}
else
{
}
if (scalar @fastq == 0)
{
	@fastq = glob("*.fastq");
}
else
{
}
if (scalar @fastq == 0)
{
	@fastq = glob("*.fastq.gz");
}
else
{
}

system("mkdir Hisat2_Readmapping_Sam");

for (my $i = 0; $i < @fastq; $i++)
{
	if (($fastq[$i] =~ /1.fq/)||($fastq[$i] =~ /1.fastq/)||($fastq[$i] =~ /_R1_/)||($fastq[$i] =~ /_1.cutadapt/)||($fastq[$i] =~ /1.non/))
	{
		if ($fastq[$i] =~ /_R2_/)
		{
			$fq2 = $fastq[$i];
			my $name = $fq2;
			$name =~ s/_2.cutadapt\S+//;
			$name =~ s/_2.fastq\S+//;
			$name =~ s/_2.fq\S+//;
			$name =~ s/^\d+[-_]//;
			print OUT ("hisat2 -p 10 --dta -x $Model -1 $fq1 -2 $fq2 -S Hisat2_Readmapping_Sam/$name.sam\n");
			system ("hisat2 -p 10 --dta -x $Model -1 $fq1 -2 $fq2 -S Hisat2_Readmapping_Sam/$name.sam");
			print OUT ("perl /var2/Lecture/Code/2.Readmapping/Hisat2_Read_Count.pl Hisat2_Readmapping_Sam/$name.sam > Hisat2_Readmapping_Sam/$name.Paired.count\n");
			system ("perl /var2/Lecture/Code/2.Readmapping/Hisat2_Read_Count.pl Hisat2_Readmapping_Sam/$name.sam > Hisat2_Readmapping_Sam/$name.Paired.count");
#			print OUT ("gzip Hisat2_Readmapping_Sam/$name.sam\n");
#			system ("gzip Hisat2_Readmapping_Sam/$name.sam");
		}
		$fq1 = $fastq[$i];
	}
	else
	{
		$fq2 = $fastq[$i];
		my $name = $fq2;
		$name =~ s/_2.cutadapt\S+//;
		$name =~ s/_2.fastq\S+//;
		$name =~ s/_2.fq\S+//;
		$name =~ s/^\d+[-_]//;
		print OUT ("hisat2 -p 10 --dta -x $Model -1 $fq1 -2 $fq2 -S Hisat2_Readmapping_Sam/$name.sam\n");
		system ("hisat2 -p 10 --dta -x $Model -1 $fq1 -2 $fq2 -S Hisat2_Readmapping_Sam/$name.sam");

		print OUT ("/var2/Lecture/Code/2.Readmapping/Hisat2_Read_Count.pl Hisat2_Readmapping_Sam/$name.sam > Hisat2_Readmapping_Sam/$name.Paired.count\n");
		system ("perl /var2/Lecture/Code/2.Readmapping/Hisat2_Read_Count.pl Hisat2_Readmapping_Sam/$name.sam > Hisat2_Readmapping_Sam/$name.Paired.count");
#		print OUT ("gzip Hisat2_Readmapping_Sam/$name.sam\n");
#		system ("gzip Hisat2_Readmapping_Sam/$name.sam");
	}
}
