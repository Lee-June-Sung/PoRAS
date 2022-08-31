use strict;

my $Fasta = $ARGV[0];

my (%bank, $id, $seq);

open (DATA, $Fasta);
while (my $stLine = <DATA>)
{
	chomp($stLine);
	my @list = split(/\s+/, $stLine);
	if ($stLine =~ />/)
	{
		$id = $list[0];
		$id =~ s/>//;
	}
	else
	{
		$bank{$id} = $bank{$id}.$stLine;
	}
}
close DATA;

my @key = sort keys %bank;

for (my $i = 0; $i < @key; $i++)
{
	my $length = length $bank{$key[$i]};
	print "$key[$i]	$length\n";
}
