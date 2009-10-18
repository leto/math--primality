#!/usr/bin/evn perl

use strict;
use warnings;
use Test::More tests => 2;

use Math::Primality;
use Math::GMPz;

### check _find_smallest_r ###
my $result = Math::Primality::_find_smallest_r(Math::GMPz->new(4099)); 
is("$result", "691", "smallest r for 4099 should be 691");
$result = Math::Primality::_find_smallest_r(Math::GMPz->new(10007)); 
is("$result", "797", "smallest r for 10007 should be 797");
