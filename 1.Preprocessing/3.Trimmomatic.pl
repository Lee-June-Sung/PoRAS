use strict;

my @file = glob("*");

my ($FQ1, $FQ2);

system("mkdir P");
system("mkdir U");


open (OUT, ">3.Trimmomatic_code.txt");
for (my $i = 0; $i < @file; $i++)
{
	if (($file[$i] =~ /1\.fastq/)||($file[$i] =~ /1\.fq/)||($file[$i] =~ /1\.gz/))
	{
		if ($file[$i] =~ /R1/)
		{
			$FQ1 = $file[$i];
		}
		elsif ($file[$i] =~ /R2/)
		{
			$FQ2 = $file[$i];
			print OUT "java -jar /home/programs/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 5 $FQ1 $FQ2 P/$FQ1\.P.fq.gz U/$FQ1\.U.fq.gz P/$FQ2\.P.fq.gz U/$FQ2\.U.fq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 -phred33\n";
			system ("java -jar /home/programs/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 5 $FQ1 $FQ2 P/$FQ1\.P.fq.gz U/$FQ1\.U.fq.gz P/$FQ2\.P.fq.gz U/$FQ2\.U.fq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 -phred33");
		}
		$FQ1 = $file[$i];
	}
	elsif (($file[$i] =~ /2\.fastq/)||($file[$i] =~ /2\.fq/)||($file[$i] =~ /2\.gz/))
	{
		$FQ2 = $file[$i];
		print OUT "java -jar /home/programs/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 5 $FQ1 $FQ2 P/$FQ1\.P.fq.gz U/$FQ1\.U.fq.gz P/$FQ2\.P.fq.gz U/$FQ2\.U.fq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 -phred33\n";
		system ("java -jar /home/programs/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 5 $FQ1 $FQ2 P/$FQ1\.P.fq.gz U/$FQ1\.U.fq.gz P/$FQ2\.P.fq.gz U/$FQ2\.U.fq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 -phred33");
	}
	
}



