#!/usr/bin/env perl

use strict;
use warnings;
use Test::More tests => 4;

use Math::Primality qw/prime_count/;
use Math::GMPz;
# http://mathworld.wolfram.com/PrimeCountingFunction.html has good test data
ok( prime_count(0) == 0, 'no primes <= 0' );
ok( prime_count(1) == 0, 'no primes <= 1' );
ok( prime_count(2) == 1, 'one prime <= 2' );
local $TODO = 'implement prime_count';
ok( prime_count(3) == 2, 'two primes <= 3' );
