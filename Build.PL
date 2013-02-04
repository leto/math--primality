use strict;
use warnings;
use Module::Build;

# NOTE: This Build.PL is to make it easier for people without dzil
# to install dependencies/etc. It may be out of sync with dist.ini
# and dist_version needs to be manually updated, for now.

my $builder = Module::Build->new(
    module_name         => 'Math::Primality',
    license             => 'perl',
    dist_author         => 'Jonathan Leto <jonathan@leto.net>',
    dist_abstract       => 'Advanced Primality Algorithms using GMP',
    dist_version        => '0.09',
    include_dirs        => '',
    build_requires => {
        'Test::More'      => 0,
        'Test::Exception' => 0,
    },
    configure_requires => { 'Module::Build' => 0.38 },
    requires => {
        'Math::GMPz'      => '0.34',
        'perl'            => '5.008',
    },
    add_to_cleanup      => [ 'Math-Primality-*' ],
    meta_merge     => {
        resources => {
            repository  => 'http://github.com/leto/math--primality/tree/master'
        },
    },
);

$builder->create_build_script();
print "Have a great day!\n";
