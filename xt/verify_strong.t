#!perl -w

use strict;
use warnings;
use lib 'lib';
use File::Spec::Functions;
use Test::More tests => 419489;
use Math::Primality qw/ is_strong_pseudoprime /;

my $bail ==<<BAIL;
You must download and unzip the file strong_pseudoprimes_less_than_1e15.txt
from http://leto.net/data/primality and put it in the xt/ directory to run
these tests. The file contains roughly half a million numbers and is 6
megabytes (2.7 gzipped)
BAIL

my $filename = 'strong_pseudoprimes_less_than_1e15.txt';

open my $fh, '<', catfile('xt',$filename) or BAIL_OUT $bail;

while(<$fh>) {
    chomp;
    ok( is_strong_pseudoprime( $_ ), "$_ is a strong pseudoprime");
}
