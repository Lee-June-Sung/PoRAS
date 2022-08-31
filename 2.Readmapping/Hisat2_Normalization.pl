use strict;

my @read_count_file = glob("Hisat2_Readmapping_Sam/*\.count");
my $gene_length = $ARGV[0];

my (%length, %last, %Read, %Read1, %last1, %Read2);

open (DATA, $gene_length);
while (my $stLine = <DATA>)
{
	chomp($stLine);
	my @list = split(/\t/, $stLine);
	$length{$list[0]} = $list[1];
}
close DATA;
my $columns = "Gene_ID";
for (my $i = 0; $i < @read_count_file; $i++)
{
	if ($read_count_file[$i] =~ /.count/)
	{
		my %bank = ();
		my @key = ();
		my $FPKM = "";
		my @list2 = split(/\//, $read_count_file[$i]);
		my $total_mapped = "";
		open (DATA, $read_count_file[$i]);
		while (my $stLine = <DATA>)
		{
			chomp($stLine);
			my @list = split(/\s+/, $stLine);
			push @key, $list[0];
			$bank{$list[0]} = $list[1];
			$total_mapped = $total_mapped + $list[1];
		}
		for (my $j = 0; $j < @key; $j++)
		{
			$Read1{$key[$j]} = $Read1{$key[$j]}."$bank{$key[$j]}\t";
			$FPKM = (1000000000*($bank{$key[$j]}/2))/($length{$key[$j]}*($total_mapped/2));	
			$last{$key[$j]} = $last{$key[$j]}."$FPKM\t";
		}
		close DATA;
	}
	my $name = $read_count_file[$i];
	$name =~ s/\.Paired\.count$//;
	$name =~ s/\.Single\.count$//;
	$name =~ s/^Hisat2_Readmapping_Sam\///;
	$columns = $columns."\t$name";
}

my @keys = sort{($a =~ /\D+\d+\D(\d+)/)[0] <=> ($b =~ /\D+\d+\D(\d+)/)[0]} keys %last;
@keys = sort{($a =~ /\D+(\d+)\D(\d+)/)[0] <=> ($b =~ /\D+(\d+)\D(\d+)/)[0]} @keys;
open (OUT, ">Total_FPKM.txt");
open (OUT3, ">Paired_Read_count.txt");
$columns =~ s/\t$//;
print OUT $columns."\n";
print OUT3 $columns."\n";
for (my $i = 0; $i < @keys; $i++)
{
	$last{$keys[$i]} =~ s/\t$//;
	$last1{$keys[$i]} =~ s/\t$//;
	$Read{$keys[$i]} =~ s/\t$//;
	$Read1{$keys[$i]} =~ s/\t$//;
	$Read2{$keys[$i]} =~ s/\t$//;
	print OUT "$keys[$i]\t$last{$keys[$i]}\n";
	print OUT3 "$keys[$i]\t$Read1{$keys[$i]}\n";
}

