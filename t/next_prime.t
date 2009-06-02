#!/usr/bin/env perl

use strict;
use warnings;
use Test::More tests => 1;

use Math::Primality qw/next_prime/;
use Math::GMPz;

### basic method handling ###
my $z = Math::GMPz->new(3);
my $p = next_prime($z);
is("$p", "5", "next_prime should handle Math::GMPz objects, 5 is the next prime after 3");
