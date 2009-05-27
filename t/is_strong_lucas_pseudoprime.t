#!perl 

use Test::More tests => 3;

use Math::Primality qw/ is_strong_lucas_pseudoprime /;
use Math::GMPz;

TODO: {
      local $TODO = "is_strong_lucas_psuedoprime is being worked on";
      my $z = Math::GMPz->new(3);
      
      ok(is_strong_lucas_pseudoprime($z), "is_strong_lucas_pseudoprime should handle Math::GMPz objects");
      ok(is_strong_lucas_pseudoprime(2), "is_strong_lucas_pseudoprime should return true for 2");
      ok(is_strong_lucas_pseudoprime(705), "is_strong_lucas_pseudoprime should return true for the first lucas pseudoprime"); 
}
