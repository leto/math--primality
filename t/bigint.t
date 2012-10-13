#!/usr/bin/env perl

use strict;
use warnings;
use Test::More;
use Math::Primality qw/:all/;
use bigint;

plan;

ok(is_prime(18446744073709551629), 'is_prime() works with bigint');
ok(is_pseudoprime(18446744073709551629), 'is_pseudoprime() works with bigint');
ok(is_pseudoprime(18446744073709551629, 3), 'is_pseudoprime(x,3) works with bigint');

ok(is_strong_lucas_pseudoprime(18446744073709551629), 'is_strong_lucas_pseudoprime() works with bigint');

done_testing;
