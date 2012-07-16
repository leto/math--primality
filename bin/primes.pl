#!/usr/bin/env perl
use strict;
use warnings;
use lib 'lib';
use Math::Primality qw/is_prime/;
$|++;

# PODNAME: primes.pl
# ABSTRACT: Print all primes between the two integers

my ($start, $end) = @ARGV;
die "USAGE:$0 start end\n" unless ($start >= 0 && $end > $start);

my $i=$start;

while ( $i++ < $end ){
    print "$i\n" if is_prime($i);
}
