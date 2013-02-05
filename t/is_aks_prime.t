#!/usr/bin/env perl

# We have to be careful with the inputs we give.  Unfortunately AKS is so
# slow that we can't have it actually prove primality using the a,n,r
# polynomial tests or we'll take a *really* long time.  With the current
# implementation including the sqrt test this means < 85991.

use strict;
use warnings;
use Test::More tests => 140;
#use Carp::Always;

use Math::GMPz qw/:mpz/;
use POSIX qw(ceil floor); # for testing _Rmpz_logbase2* functions

BEGIN {
  use_ok ('Math::Primality::AKS' );
}
use Math::Primality::AKS qw/is_aks_prime/;

my $z = Math::GMPz->new(3);
ok( is_aks_prime($z), "is_aks_prime should handle Math::GMPz objects, three is prime" );
ok( is_aks_prime(2), '2 is prime');
ok(!is_aks_prime(1), '1 is not prime');
ok(!is_aks_prime(0), '0 is not prime');
ok(!is_aks_prime(-1), '-1 is not prime');

# What should it do for x < 0 ?
#ok(!is_aks_prime(-2), '-2 is not prime');
ok( !is_aks_prime(20), "20 is not prime");

# powers of 2 are never prime
for my $k (1..20) {
        Rmpz_ui_pow_ui($z, 2, ++$k );
        ok(!is_aks_prime($z), "2**$k=$z is not prime");
}

my @small_primes = qw/
  5   7  11  13  17  19  23  29  31  37  41  43  47  53  59  61  67  71  73
 79  83  89  97 101 103 107 109 113 127 131 137 139 149 151 157 163 167 173
179 181 191 193 197 199 211 223 227 229 233 239 241 251 257 263 269 271 277
281 283 293 307 311 313 317 331 337 
/;

map { ok(is_aks_prime($_), "$_ is an AKS prime") } @small_primes;

my @carmichael = qw/561 1105 1729 2465 2821 6601 8911
                    10585 15841 29341 41041 46657 52633
                    62745 63973 75361 101101
                    999838193331601
                    999840927672001
                    999851057445241
                    999878556600001
                    999885684921481
                    999895175363161
                    999902676805201
                    999907821232321
                    999919121100481
                    999922265173441
/;
map { ok(!is_aks_prime($_), "Carmichael Number $_ is not an AKS prime") } @carmichael;

# First 20 psp(2)'s
map { ok(!is_aks_prime($_), "Pseudoprime (base 2) $_ is not an AKS prime" ) } qw/
 341 561 645 1105 1387 1729 1905 2047
 2465 2701 2821 3277 4033 4369 4371
 4681 5461 6601 7957 8321
/;

