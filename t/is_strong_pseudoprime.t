#!perl -T

use Test::More tests => 58;

use Math::Primality qw/ is_strong_pseudoprime /;

ok(!is_strong_pseudoprime(-1),'-1 is not a spsp' );
ok(!is_strong_pseudoprime(0),'0 is not a spsp' );
ok(!is_strong_pseudoprime(1),'1 is not a spsp' );
ok( is_strong_pseudoprime(2,3),'2 is a spsp(3)' );
ok( is_strong_pseudoprime(3),'3 is a spsp'  );
ok(!is_strong_pseudoprime(4),'4 is not a spsp' );
ok( is_strong_pseudoprime( 1093**2 ), '1093**2 is a spsp');
ok(!is_strong_pseudoprime(1000), '1000 is not a spsp');


while(<DATA>) {
    chomp;
    ok( is_strong_pseudoprime( $_ ), "$_ is a spsp");
}


__DATA__
2047
3277
4033
4681
8321
15841
29341
42799
49141
52633
65281
74665
80581
85489
88357
90751
104653
130561
196093
220729
233017
252601
253241
256999
271951
280601
314821
357761
390937
458989
476971
486737
489997
514447
580337
635401
647089
741751
800605
818201
838861
873181
877099
916327
976873
983401
1004653
1016801
1023121
1082401
