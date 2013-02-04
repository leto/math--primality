#!/usr/bin/env perl
use strict;
use warnings;
use lib 'lib';
use Math::Prime::Util;
use Math::Prime::Util::GMP;
use Math::Primality;
use Benchmark qw/:all/;
my $count = shift || -2;
srand(29);  # So we have repeatable results

test_at_digits($_, 1000) for (5, 15, 25, 50, 200);

sub test_at_digits {
  my($digits, $numbers) = @_;
  die "Digits must be > 0" unless $digits > 0;

  my @nums = map { Math::Prime::Util::random_ndigit_prime($digits)+2 } 1 .. $numbers;
  print "is_prime for $numbers random $digits-digit numbers\n";

  cmpthese($count,{
    'MP'   => sub { Math::Primality::is_prime($_) for @nums; },
    'MPU'  => sub { Math::Prime::Util::is_prob_prime($_) for @nums; },
  });
}
