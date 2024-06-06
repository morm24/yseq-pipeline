use strict;
use warnings;
use Getopt::Long;

GetOptions(
	'column=i' => \(my $column = 0),
	'filename=s' => \(my $filename = "")
) or die "Invalid options passed to $0\n";

my $sum = 0;
#print "\n" . $column . "\n";
#print $filename . "\n";

open(my $fh, '<', $filename);

while( my $line = <$fh>) {
	chomp($line);
	my @spl = split('\t',$line);
	$sum += $spl[$column];
}
close $fh;

# open my $fw, '>&', STDOUT or die "bleah: $!";
# print $fw $sum;
print "$sum";
exit(1);
