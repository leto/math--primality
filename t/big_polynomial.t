#!/usr/bin/evn perl

use strict;
use warnings;
use Test::More;
#use Carp::Always;

use Math::GMPz qw/:mpz/;

BEGIN {
  use_ok ('Math::Primality::BigPolynomial' );
}

sub new_poly {
 return Math::Primality::BigPolynomial->new([ map { Math::GMPz->new($_) } @_ ]);
}

my $b = Math::Primality::BigPolynomial->new([1,3,7]);
isa_ok($b, 'Math::Primality::BigPolynomial');

is($b->getCoef(0), 1, 'coef(0) is 1');
is($b->getCoef(1), 3, 'coef(1) is 3');
is($b->getCoef(2), 7, 'coef(2) is 7');
# this is a bit wonky
cmp_ok($b->getCoef(3),'==', Math::GMPz->new(0), 'coef(3) is 0');
is($b->getCoef(-1), undef, 'coef(-1) is undef');
is($b->getCoef(-10), undef, 'coef(-10) is undef');

# Now do some testing for multiplication
{
  my $p1 = new_poly(3,1);
  my $p2 = new_poly(2,1);
  is($p1->copy()->mulmod($p2,1024,1024)->string,
     "x^2 + 5x + 6",
     "(x+3) * (x+2) = x^2 + 5x + 6");
}
{
  my $p1 = new_poly(4, 0, 2, 1);
  my $p2 = new_poly(1, 1, 0, 2);
  is($p1->copy()->mulmod($p2,1024,1024)->string,
     "2x^6 + 4x^5 + x^4 + 11x^3 + 2x^2 + 4x + 4",
     "(x^3 + 2x^2 + 4) * (2x^3 + x + 1)");
  is($p1->copy()->powmod(3, 1024, 1024)->string,
     "x^9 + 6x^8 + 12x^7 + 20x^6 + 48x^5 + 48x^4 + 48x^3 + 96x^2 + 64",
     "(x^3 + 2x^2 + 4) ^ 3");
}
{
  # http://reference.wolfram.com/mathematica/tutorial/PolynomialsModuloPrimes.html
  my $p = new_poly(1,1);
  is($p->copy->powmod(6, 1024, 1024)->string,
     "x^6 + 6x^5 + 15x^4 + 20x^3 + 15x^2 + 6x + 1",
     "(1+x)^6");
  is($p->copy->powmod(6, 2, 1024)->string,
     "x^6 + x^4 + x^2 + 1",
     "(1+x)^6 modulo 2");
}

done_testing;

