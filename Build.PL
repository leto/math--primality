use strict;
use warnings;
use Module::Build;

my $builder = Module::Build->new(
    module_name         => 'Math::Primality',
    license             => 'perl',
    dist_author         => 'Jonathan Leto <jonathan@leto.net>',
    dist_abstract       => 'Advanced Primality Algorithms using GMP',
    dist_version_from   => 'lib/Math/Primality.pm',
    include_dirs        => '',
    build_requires => {
        'Test::More' => 0,
        'Test::Exception' => 0,
    },
    requires => {
        'Math::GMPz' => 0.26,
    },
    add_to_cleanup      => [ 'Math-Primality-*' ],
    create_makefile_pl  => 'traditional',
);

$builder->create_build_script();
print "Have a great day!\n";
