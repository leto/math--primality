#!perl -w

use strict;
use warnings;
use lib 'lib';
use File::Spec::Functions;
use Test::More tests => 419489;
use Math::Primality qw/ is_strong_pseudoprime /;

my $bail = <<BAIL;
You must download and unzip the file strong_pseudoprimes_less_than_1e15.txt
from http://leto.net/data/primality and put it in the xt/ directory to run
these tests. The file contains roughly half a million numbers and is 6
megabytes (2.7 gzipped)

Example commands to get this data:

cd xt/
wget http://leto.net/data/primality/spsp_base2_1e15.txt.gz
gunzip spsp_base2_1e15.txt.gz

BAIL

my $filename = catfile(qw/xt spsp_base2_1e15.txt/);

open my $fh, '<', $filename or bail();

sub bail
{
    diag $bail;
    exit(1);
}

while(<$fh>) {
    chomp;
    ok( is_strong_pseudoprime( $_ ), "$_ is a strong pseudoprime");
}

close $fh;
