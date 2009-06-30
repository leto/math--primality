#!/usr/bin/env perl

use strict;
use warnings;
use lib 'lib';
use File::Spec::Functions;
#use Test::More tests => 1270607;
print "1..1270607\n";
use Math::Primality qw/ is_prime /;
$|++; #flush the output buffer after every write() or print() function

my $bail = <<BAIL;

RUH-ROH, Shaggy.  You need primes_less_than_2e10.txt. 

BAIL

my $filename = catfile(qw/xt primes_less_than_2e10.txt/);

open my $fh, '<', $filename or bail();

sub bail
{
    diag $bail;
    exit(1);
}

my $t = 1;
while(<$fh>) {
    chomp;
    $_ = trim($_);
    print is_prime( $_ ) ? "ok $t - $_ is a prime\n" : "not ok $t - $_ is a prime\n";
    $t++;
}

sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}

close $fh;
