#!/usr/bin/env perl

use strict;
use warnings;
use Test::More tests => 12;

use Math::Primality qw/ is_strong_lucas_pseudoprime/;
use Math::GMPz;

ok(is_strong_lucas_pseudoprime(2), "is_strong_lucas_pseudoprime should return true for 2");

TODO: {
      local $TODO = "is_strong_lucas_psuedoprime is being worked on";
      my $z = Math::GMPz->new(3);

      ok(is_strong_lucas_pseudoprime($z), "is_strong_lucas_pseudoprime should handle Math::GMPz objects");
      ok(is_strong_lucas_pseudoprime(705), "is_strong_lucas_pseudoprime should return true for the first lucas pseudoprime"); 
};

ok(!is_strong_lucas_pseudoprime(9), 'is_strong_lucas_pseudoprime deals with perfect squares');

### Test _find_dpq_selfridge ###
ok(
    eq_array( [Math::Primality::_find_dpq_selfridge(1993) ], [ 5, 1, -1 ]),
    "_find_dpq_selfridge should return (5, 1, -1) for 1993"
);

ok(
    eq_array( [Math::Primality::_find_dpq_selfridge(1759) ], [-11, 1, 3]),
    "_find_dpq_selfridge should return (-11, 1, 3) for 1759"
);

### Test _find_s_d ###
{
  my ($s, $d) = Math::Primality::_find_s_d(7);
  ok ($s == 1, "_find_s_d should return 7 = 3 * 2^1 + 1");
  is ("$d", "3",  "_find_s_d should return 7 = 3 * 2^1 + 1");  
}

{
  my ($s, $d) = Math::Primality::_find_s_d(17);
  ok ($s == 4, "_find_s_d should return 17 = 1 * 2^4 + 1");
  is ("$d", "1", "_find_s_d should return 17 = 1 * 2^4 + 1");
}

{
  my ($s, $d) = Math::Primality::_find_s_d(53525);
  ok ($s == 2, "_find_s_d should return 53525 = 13381 * 2^2 + 1");
  is ("$d", "13381", "_find_s_d should return 53525 = 13381 * 2^2 + 1");
}
