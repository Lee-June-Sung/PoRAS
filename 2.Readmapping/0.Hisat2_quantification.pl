use strict;

my $DB = $ARGV[0];
my $com;

open (OUT, ">0.Commend.out");

system ("hisat2-build $DB $DB");
system ("perl /var2/Lecture/Code/2.Readmapping/Hisat2_Read_Mapping.pl $DB");
system ("perl /var2/Lecture/Code/2.Readmapping/Gene_Length.pl $DB > $DB.length");
system ("perl /var2/Lecture/Code/2.Readmapping/Hisat2_Normalization.pl $DB.length");
system ("pypy /var2/Lecture/Code/2.Readmapping/sam_to_mapped_read.py Hisat2_Readmapping_Sam");

print OUT "hisat2-build $DB $DB\n";
print OUT "perl /var2/Lecture/Code/2.Readmapping/Hisat2_Read_Mapping.pl $DB\n";
print OUT "perl /var2/Lecture/Code/2.Readmapping/Gene_Length.pl $DB > $DB.length\n";
print OUT "perl /var2/Lecture/Code/2.Readmapping/Hisat2_Normalization.pl $DB.length\n";
print OUT "pypy /var2/Lecture/Code/2.Readmapping/sam_to_mapped_read.py Hisat2_Readmapping_Sam\n";



