use strict;

my @file = glob("*");

my ($FQ1, $FQ2);
open (OUT, ">2.Cut_Adapt.txt");
system("mkdir Cutadapt");
for (my $i = 0; $i < @file; $i++)
{
	if (($file[$i] =~ /1\.fastq/)||($file[$i] =~ /1\.fq/))
	{
		if ($file[$i] =~ /R1/)
		{
			$FQ1 = $file[$i];
		}
		elsif ($file[$i] =~ /R2/)
		{
			$FQ2 = $file[$i];
			print OUT "cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o Cutadapt/$FQ1.NoAdapt.fq.gz -p Cutadapt/$FQ2.NoAdapt.fq.gz $FQ1 $FQ2\n";
			system("cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o Cutadapt/$FQ1.NoAdapt.fq.gz -p Cutadapt/$FQ2.NoAdapt.fq.gz $FQ1 $FQ2");	
		}
		$FQ1 = $file[$i];
	}
	elsif (($file[$i] =~ /2\.fastq/)||($file[$i] =~ /2\.fq/))
	{
		$FQ2 = $file[$i];
		print OUT "cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o Cutadapt/$FQ1.NoAdapt.fq.gz -p Cutadapt/$FQ2.NoAdapt.fq.gz $FQ1 $FQ2\n";
		system("cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o Cutadapt/$FQ1.NoAdapt.fq.gz -p Cutadapt/$FQ2.NoAdapt.fq.gz $FQ1 $FQ2");
	}
}
