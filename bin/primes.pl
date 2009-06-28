#!/usr/bin/env perl
use strict;
use warnings;
use lib 'lib';
use Math::Primality qw/is_prime/;
$|++;

my ($start, $end) = @ARGV;
die "USAGE:$0 start end\n" unless ($start >= 0 && $end > $start);

my $i=$start;

while ( $i++ <= $end ){
    print "$i\n" if is_prime($i);
}
