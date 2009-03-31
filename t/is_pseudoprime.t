#!perl -T

use Test::More tests => 48;

use Math::Primality qw/ is_pseudoprime /;


ok( is_pseudoprime(7,2) , '7 is a base-2 psp');
ok( is_pseudoprime(7) , 'is_psp defaults to base 2');
ok(!is_pseudoprime(2), '2 is not a psp' );
ok(!is_pseudoprime(4), '4 is not a psp' );

# These are the 10 largest Carmichael numbers <= 1e15
my @carmichael = qw/
                    999838193331601
                    999840927672001
                    999851057445241
                    999878556600001
                    999885684921481
                    999895175363161
                    999902676805201
                    999907821232321
                    999919121100481
                    999922265173441
                    /;
# non of the above numbers have these as factors
my @bases = (2, 3, 5, 23);
for my $base (@bases) {
    map {
        ok( is_pseudoprime($_,$base), "$_ is a Carmicheal number and hence a psp($base)");
    } @carmichael;
};

ok(!is_pseudoprime( 999895175363161, 7 ), ' 7 is a factor of 999895175363161 and hence it is not a psp(7)');

map {
ok(!is_pseudoprime( $_ , 13 ), "13 is a factor of $_ and hence it is not a psp(13)") } (999838193331601,999878556600001,999902676805201) ;
