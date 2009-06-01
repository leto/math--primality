#!/usr/bin/env perl

use strict;
use warnings;
use Test::More tests => 5;

use Math::Primality qw/ is_strong_lucas_pseudoprime _find_dpq_selfridge /;
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
    eq_array( [_find_dpq_selfridge(1993) ], [ 5, 1, -1 ]),
    "_find_dpq_selfridge should return (5, 1, -1) for 1993"
);
