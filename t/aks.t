#!/usr/bin/evn perl

use strict;
use warnings;
use Test::More tests => 3;

use Math::Primality;
use Math::GMPz;
use POSIX qw(ceil floor); # for testing _Rmpz_logbase2* functions

### check _find_smallest_r ###
my $result = Math::Primality::_find_smallest_r(Math::GMPz->new(4099)); 
is("$result", "691", "smallest r for 4099 should be 691");
$result = Math::Primality::_find_smallest_r(Math::GMPz->new(10007)); 
is("$result", "797", "smallest r for 10007 should be 797");
### log base 2 for testing my _Rmpz_logbase2* functions
sub log2 {
  my $n = shift;
  return log($n)/log(2);
}
### check _Rmpz_logbase2cl ###
$a = ceil(log2(100));
$b = Math::Primality::_Rmpz_logbase2cl(Math::GMPz->new(100));
is("$a", "$b", "log2(100) == _Rmpz_logbase2cl(100)");
### check _Rmpz_logbase2fl ###
$a = floor(log2(100));
$b = Math::Primality::_Rmpz_logbase2cl(Math::GMPz->new(100));
is("$a", "$b", "log2(100) == _Rmpz_logbase2cl(100)");
