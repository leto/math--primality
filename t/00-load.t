#!/usr/bin/env perl
use strict;

use Test::More tests => 1;

BEGIN {
    use_ok( 'Math::Primality' );
}

diag( "Testing Math::Primality $Math::Primality::VERSION, Perl $], $^X" );
