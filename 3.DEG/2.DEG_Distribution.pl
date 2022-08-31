use strict;

my $file = $ARGV[0];
my $out = $ARGV[1];
my $FC = $ARGV[2];

open (OUT, ">$out.out");
open (OUT2, ">$out.table");
my @files = split(/\,/, $file);
my $cnt = scalar @files;

my (%check1, %bank, %Gene);
print OUT $_."\t" foreach (@files);
print OUT "\n";

for (my $i = 0; $i < $cnt; $i++)
{
	open(DATA, $files[$i]);
	while (my $stLine = <DATA>)
	{
		chomp($stLine);
		$stLine =~ s/"//g;
		my @list = split(/\,/, $stLine);
		$Gene{$list[1]}++;
	}
}
my @geneID = sort keys %Gene;

for (my $j = 0; $j < $cnt; $j++)
{
	my %ch = ();
	open (DATA, $files[$j]);
	while (my $stLine = <DATA>)
	{
		chomp($stLine);
		$stLine =~ s/"//g;
#		$stLine =~ s/NA/0/g;
		my @list = split(/\,/, $stLine);
#		if (($list[8] > 20)||($list[9] > 20)||($list[10] > 20)||($list[11] > 20)||($list[12] > 20)||($list[13] > 20))
#		{
		if ($list[7] ne "NA")
		{
			if ($list[7] < 0.05)
			{
				#if ($list[6] < 0.05)
				#{
					if ($list[3] >= $FC)
					{	
						$check1{$list[1]} = $check1{$list[1]}."+$j ";
#						$bank{$list[1]} = $bank{$list[1]}."$list[3]	$list[6]	$list[7]	";
						$ch{$list[1]} = "OK";
					}
					elsif ($list[3] <= -$FC)
					{
						$check1{$list[1]} = $check1{$list[1]}."-$j ";
#						$bank{$list[1]} = $bank{$list[1]}."$list[3]	$list[6]	$list[7]	";
						$ch{$list[1]} = "OK";
					}
				#}
			}
		}	
#		}
		$bank{$list[1]} = $bank{$list[1]}."$list[3]	$list[6]	$list[7]	";
		
	}
	close DATA;

	for (my $i = 0; $i < @geneID; $i++)
	{
		if ($ch{$geneID[$i]} eq "")
		{
#			$bank{$geneID[$i]} = $bank{$geneID[$i]}."0\t0\t0\t";
		}
	}
print OUT2 "\t$files[$j]\_Log2FC\t$files[$j]\_P-value\t$files[$j]\_adj.P-value";
}
my @key = sort keys %check1;

#open (OUT, ">Check.out");
print OUT2 "\tGroup_Information\n";
for (my $i = 0; $i < @key; $i++)
{
	print OUT "$key[$i]	$check1{$key[$i]}\n";
	my @list = split(/\s+/, $bank{$key[$i]});
	my $check = "";
	print OUT2 "$key[$i]	$bank{$key[$i]}$check1{$key[$i]}\n";
}

