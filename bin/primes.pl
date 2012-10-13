#!/usr/bin/env perl
use strict;
use warnings;
use lib 'lib';
use Math::Primality qw/is_prime next_prime/;
use Math::GMPz;
$|++;

# PODNAME: primes.pl
# ABSTRACT: Print all primes between the two integers

my ($start, $end) = @ARGV;
die "USAGE:$0 start end\n" unless (@ARGV == 2);

die "$start isn't a positive integer" if $start =~ tr/0123456789//c;
die "$end isn't a positive integer" if $end =~ tr/0123456789//c;

$start = Math::GMPz->new("$start");
$end   = Math::GMPz->new("$end");
$start = next_prime($start) unless is_prime($start);
while ($start <= $end) {
    print "$start\n";
    $start = next_prime($start);
}
