use strict;
use warnings;
use Getopt::Long;

GetOptions(
	'multiplicands=s' => \(my $multiplicandsString = (0)),
	'round' => \(my $round = 0)
) or die "Invalid options passed to $0\n";

my $product = 1;

my @multiplicands = split(',', $multiplicandsString);

for my $multiplicand (@multiplicands) {
	$product *= $multiplicand;
}
if ($round) {
	$product = int($product + 0.5);
}
print $product;
