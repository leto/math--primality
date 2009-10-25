#!/usr/bin/evn perl

use strict;
use warnings;
use Test::More tests => 14;

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
for (my $i = 1; $i < 7; $i++) {
  my $n = $i * 100;
  my $a = ceil(log2($n));
  my $b = Math::Primality::_Rmpz_logbase2cl(Math::GMPz->new($n));
  is("$a", "$b", "log2($n) == _Rmpz_logbase2cl($n)");
}
### check _Rmpz_logbase2fl ###
for (my $i = 1; $i < 7; $i++) {
  my $n = $i * 100;
  $a = floor(log2($n));
  $b = Math::Primality::_Rmpz_logbase2fl(Math::GMPz->new($n));
  is("$a", "$b", "log2($n) == _Rmpz_logbase2fl($n)");
}
