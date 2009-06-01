package Math::Primality;
use warnings;
use strict;
use Data::Dumper;
use Math::GMPz qw/:mpz/;
use base 'Exporter';
our $DEBUG = 0;

use constant GMP => 'Math::GMPz';

=head1 NAME

Math::Primality - Advanced Primality Algorithms using GMP

=head1 VERSION

Version 0.02

=cut

our $VERSION = '0.03_01';

our @EXPORT_OK = qw/is_pseudoprime is_strong_pseudoprime is_prime is_strong_lucas_pseudoprime _find_dpq_selfridge /;

our %EXPORT_TAGS = ( all => \@EXPORT_OK );

=head1 SYNOPSIS

    use Math::Primality qw/:all/;

    my $t1 = is_pseudoprime($x,$base);
    my $t2 = is_strong_pseudoprime($x);

=head1 EXPORT

=head1 FUNCTIONS

=head2 is_pseudoprime($n,$b)

Returns true if $n is a base $b pseudoprime, otherwise false.  The variable $n
should be a Perl integer or Math::GMPz object.

The default base of 2 is used if no base is given. Base 2 pseudoprimes are often called Fermat pseudoprimes.

    if ( is_pseudoprime($n,$b) ) {
        # it's a pseudoprime
    } else {
        # not a psuedoprime
    }

=cut

sub debug {
    warn $_[0] if $ENV{DEBUG} or $DEBUG;
}

sub is_pseudoprime($;$)
{
    my ($n, $base) = @_;
    return 0 unless $n;
    $base ||= 2;
    # we should check if we are passed a GMPz object
    $base   = GMP->new($base);
    $n      = GMP->new($n);

    # if $n and $base are not coprime, than $base is a factor of $n
    # $base > 2 && ( Math::BigInt::bgcd($n,$base) != 1 ) && return 0;

    my $m   = _copy($n);
    Rmpz_sub_ui($m, $m, 1);              # m = n - 1

    my $mod = _copy($base);
    Rmpz_powm($mod, $base, $m, $n );
    return ! Rmpz_cmp_ui($mod, 1);       # pseudoprime if $mod = 1
}

sub _copy($)
{
    my ($n) = @_;
    return GMP->new($n);
}

=head2 is_strong_pseudoprime($n,$b)

Returns true if $n is a base $b strong pseudoprime, false otherwise.  The variable $n should be a Perl integer
or a Math::GMPz object. Strong psuedoprimes are often called Miller-Rabin pseudoprimes.

The default base of 2 is used if no base is given.

    if ( is_strong_pseudoprime($n,$b) ) {
        # it's a strong pseudoprime
    } else {
        # not a strong psuedoprime
    }

=cut

sub is_strong_pseudoprime($;$)
{
    my ($n, $base) = @_;

    $base ||= 2;
    # we should check if we are passed a GMPz object
    $base   = GMP->new($base);
    $n      = GMP->new($n);

    my $cmp = Rmpz_cmp_ui($n, 2 );
    return 1 if $cmp == 0;
    return 0 if $cmp < 0;

    # unnecessary but faster if $n is even
    return 0 if Rmpz_even_p($n);

    my $m   = _copy($n);
    Rmpz_sub_ui($m,$m,1);

    my $s   = Rmpz_scan1($m,0);
    my $d   = GMP->new(0);

    Rmpz_tdiv_q_2exp($d, $m,$s);
    debug "m=$m, s=$s, d=$d";

    my $residue = GMP->new(0);
    Rmpz_powm($residue, $base,$d, $n);
    debug "$base^$d % $n = $residue";

    # if $base^$d = +-1 (mod $n) , $n is a strong pseudoprime

    if ( Rmpz_cmp_ui( $residue,1) == 0 ) {
        debug "found $n as spsp since $base^$d % $n == $residue == 1\n";
        return 1;
    }

    if ( Rmpz_cmp($residue,$m) == 0 ) {
        debug "found $n as spsp since $base^$d % $n == $residue == $m\n";
        return 1;
    }

    map {
        # successively square $residue, $n is a strong psuedoprime
        # if any of these are congruent to -1 (mod $n)
        Rmpz_mul($residue,$residue,$residue);
        debug "$_: r=$residue";

        my $mod = _copy($residue);
        Rmpz_mod($mod, $mod,$n);
        debug "$_:$residue % $n = $mod ";
        $mod = Rmpz_cmp($mod, $m);

        if ($mod == 0) {
            debug "$_:$mod == $m => spsp!";
            return 1;
        }
    } ( 1 .. $s-1 );

    return 0;
}

sub is_strong_lucas_pseudoprime($)
{
    my ($n) = @_;
    $n      = GMP->new($n);
    # weed out all perfect squares
    if ( Rmpz_perfect_square_p($n) ) {
        return 0;
    }
    # we also need to weed out all N < 3 and all even N 
    my $cmp = Rmpz_cmp_ui($n, 2 );
    return 1 if $cmp == 0;
    return 0 if $cmp < 0;
    return 0 if Rmpz_even_p($n);
    # determine Selfridge parameters D, P and Q
    my ($d, $p, $q) = _find_dpq_selfridge($n);
    # translate code which computes lucas numbers in iStrongLucasSelfridge
    return 0;
}

#selfridge's method for finding the tuple (D,P,Q) for is_strong_lucas_pseudoprime
sub _find_dpq_selfridge($) {
  my $n = GMP->new($_[0]);
  my ($d,$sign,$wd) = (5,1,0);
  my $gcd = Math::GMPz->new;

  # determine D
  while (1) {
    $wd = $d * $sign;

    Rmpz_gcd_ui($gcd, $n, abs $wd);
    if ($gcd > 1 && Rmpz_cmp_ui($n, $gcd) > 0) {
      debug "1 < $gcd < $n => $n is composite with factor $wd";
      return 0;
    }
    my $j = Rmpz_jacobi(GMP->new($wd), $n);
    if ($j == -1) {
      debug "Rmpz_jacobi($wd, $n) == -1 => found D";
      last; 
    }
    # didn't find D, increment and swap sign
    $d += 2;
    $sign = -$sign;
    ### TODO ###
    # need code to make sure we don't overflow $d 
    ### TODO ###
  }
  # P = 1
  my ($p,$q) = (1,0);
  {
      use integer;
      # Q = (1 - D) / 4
      $q = (1 - $wd) / 4;
  }
  return ($wd, $p, $q);
}

#alternate method for finding the tuple (D,P,Q) for is_strong_lucas_pseudoprime
sub _find_dpq_alternate($) {

}

=head1 AUTHOR

Jonathan Leto, C<< <jonathan at leto.net> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-math-primality at
rt.cpan.org>, or through the web interface at
L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Math::Primality>.  I will be
notified, and then you'll automatically be notified of progress on your bug as I
make changes.


=head1 THANKS

The algorithms in this module have been ported from the C source code in
bpsw1.zip by Thomas R. Nicely, available at http://www.trnicely.net/misc/bpsw.html
or in the spec/bpsw directory of the Math::Primality source code. Without his
research this module would not exist.


=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Math::Primality


You can also look for information at:

=over 4

=item * Math::Primality on Github

L<http://github.com/leto/math--primality/tree/master>

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

exp(0); # End of Math::Primality
