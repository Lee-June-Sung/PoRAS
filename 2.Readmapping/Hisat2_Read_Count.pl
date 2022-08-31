use strict;

my $sam = $ARGV[0];
	my $name = $sam;
	$name =~ s/\.sam$//;
	my (%bank, @key);
	open (DATA, $sam);
	while (my $stLine = <DATA>)
	{
		chomp($stLine);
		my @list = split(/\s+/, $stLine);
		if ($list[0] =~ /^@/)
		{
	#		print $list[1]."\n";
			if($list[1] =~ /SN:/)
			{
				$list[1] =~ s/SN://;
	#			print $list[1]."\n";
				$bank{$list[1]} = 0;
				push @key, $list[1];
			}
		}
		else
		{
			$bank{$list[2]}++;
			if (($list[6] ne "\=")||($list[6] ne "\*"))
			{
				$bank{$list[6]}++;
			}
		}
	}
	close DATA;
	
	for (my $i = 0; $i < @key; $i++)
	{
		print "$key[$i]	$bank{$key[$i]}\n";
	}
