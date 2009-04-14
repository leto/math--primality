#!perl -T

use Test::More tests => 9;

use Math::Primality qw/ is_strong_pseudoprime /;

ok(!is_strong_pseudoprime(-1),'-1 is not a spsp' );
ok(!is_strong_pseudoprime(0),'0 is not a spsp' );
ok(!is_strong_pseudoprime(1),'1 is not a spsp' );
ok( is_strong_pseudoprime(2,3),'2 is a spsp(3)' );
ok( is_strong_pseudoprime(3),'3 is a spsp'  );
ok(!is_strong_pseudoprime(4),'4 is not a spsp' );

ok( is_strong_pseudoprime( 2047 ), '2047 is a spsp');
ok( is_strong_pseudoprime( 1093**2 ), '1093**2 is a spsp');
ok( is_strong_pseudoprime( 561, 3 ), '561 is a spsp(3)');
