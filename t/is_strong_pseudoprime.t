#!perl -T

use Test::More tests => 4;

use Math::Primality qw/ is_strong_pseudoprime /;

ok(!is_strong_pseudoprime(0) );
ok(!is_strong_pseudoprime(1) );
ok( is_strong_pseudoprime(2) );
ok( is_strong_pseudoprime(3) );

