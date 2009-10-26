#!/usr/bin/evn perl

use strict;
use warnings;
use Test::More tests => 16;

use Math::Primality;
use Math::GMPz;
use POSIX qw(ceil floor); # for testing _Rmpz_logbase2* functions

### check _find_smallest_r ###
my $smallest_r = Math::Primality::_find_smallest_r(Math::GMPz->new(4099)); 
is("$smallest_r", "691", "smallest r for 4099 should be 691");
$smallest_r = Math::Primality::_find_smallest_r(Math::GMPz->new(10007)); 
is("$smallest_r", "797", "smallest r for 10007 should be 797");
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
### check _calculate_max_a ###
my $max_a = Math::Primality::_calculate_max_a(Math::GMPz->new(797), Math::GMPz->new(10007));
is("$max_a", "789", "max a of 10007 (with r of 797) should be 789");
$max_a = Math::Primality::_calculate_max_a(Math::GMPz->new(691), Math::GMPz->new(4099));
is("$max_a", "682", "max a of 4099 (with r of 691) should be 682");
