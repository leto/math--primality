#!perl -w
use strict;
use warnings;
use lib 'lib';
use Math::Primality qw/is_strong_pseudoprime/;
$|++;

my ($base, $start, $end) = @ARGV;
die "USAGE:$0 base start end\n" unless ($base && $start >= 0 && $end > $start);

my $i=$start;

# This currently includes primes as well, which need to be sifted out with 
# !is_prime when it is written
print "Generating spsp($base)\n";
while ( $i++ <= $end ){
    print "$i\n" if is_strong_pseudoprime($i,$base); # && !is_prime($i)
}



