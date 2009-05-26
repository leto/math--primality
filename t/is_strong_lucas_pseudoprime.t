#!perl -T

use Test::More tests => 1;

use Math::Primality qw/ is_strong_lucas_pseudoprime /;

ok(is_strong_lucas_pseudoprime(2) );
