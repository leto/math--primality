#!/usr/bin/env perl

use strict;
use warnings;
use Test::More tests => 3;

use Math::Primality qw/is_prime/;
use Math::GMPz;

### basic method handling ###
my $z = Math::GMPz->new(3);
ok( is_prime(3), "is_prime should handle Math::GMPz objects, three is prime" );
ok( is_prime(2), "is_prime should handle 2 as a prime");
ok( !is_prime(20), "is_prime should even numbers");

TODO: {
      local $TODO = "is_prime is being worked on";

      ### test small numbers (<1000) ###
      ### test known Miller-Rabin psuedoprimes base 2 ###
      ### test Carmichael numbers ###
      ### test Lucas psuedoprimes ###
};
