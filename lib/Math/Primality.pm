package Math::Primality;

use Math::BigInt qw/bgcd/;
use base 'Exporter';
use warnings;
use strict;

=head1 NAME

Math::Primality - Various Primality Algorithms

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';

our @EXPORT_OK = qw/is_pseudoprime/;

our %EXPORT_TAGS = ( all => \@EXPORT_OK );

=head1 SYNOPSIS

    use Math::Primality;
    use Math::BigInt;

    my $t1 = is_pseudoprime($x,$base);
    my $t2 = is_lucas_pseudoprime($x);

=head1 EXPORT

=head1 FUNCTIONS

=head2 is_pseudoprime($n,$b)

Returns true if $n is a base $b pseudoprime, otherwise false.
The default base of 2 is used if no base is given. Base 2 pseudoprimes are often called Fermat pseudoprimes.

    if ( is_pseudoprime($n,$b) ) {
        ...
    } else {
        ...
    }

=cut

sub is_pseudoprime
{
    my ($n, $base) = @_;
    # force to BigInts for now
    $base ||= 2;
    $base   = Math::BigInt->new("$base");
    $n      = Math::BigInt->new("$n");

    # if $n and $base are not coprime, than $base is a factor of $n
    # $base > 2 && ( Math::BigInt::bgcd($n,$base) != 1 ) && return 0;

    my $m   = $n->copy->bdec;            # m = n - 1
    my $mod = $base->copy;
    $mod->bmodpow($m,$n);     # (base**exp) (mod n)
    return $mod == 1;
}

=head2 is_strong_pseudoprime($n,$b)

Returns true if $n is a base $b strong pseudoprime, false otherwise.

=cut

sub is_strong_pseudoprime
{
    my ($n, $base) = @_;
    # force to BigInts for now
    $base ||= 2;
    $base   = Math::BigInt->new("$base");
    $n      = Math::BigInt->new("$n");
    my $m   = $n->copy->bdec();

    # need access to mpz_scan1

    return 0;
}
=head1 AUTHOR

Jonathan Leto, C<< <jonathan at leto.net> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-math-primality at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Math::Primality>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Math::Primality


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=Math::Primality>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/Math::Primality>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/Math::Primality>

=item * Search CPAN

L<http://search.cpan.org/dist/Math::Primality>

=back


=head1 ACKNOWLEDGEMENTS


=head1 COPYRIGHT & LICENSE

Copyright 2009 Jonathan Leto, all rights reserved.

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.


=cut

1; # End of Math::Primality
