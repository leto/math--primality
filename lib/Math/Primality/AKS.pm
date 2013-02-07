package Math::Primality::AKS;
use warnings;
use strict;

use Math::GMPz qw/:mpz/;
use Math::Primality::BigPolynomial;

use base 'Exporter';
use Carp qw/croak/;

use constant DEBUG => 0;

use constant GMP => 'Math::GMPz';

# ABSTRACT: Check for primality using the AKS (Agrawal-Kayal-Saxena) algorithm

=head1 NAME

Math::Primality::AKS - Check for primes with AKS

=cut

our @EXPORT_OK = qw/is_aks_prime/;

our %EXPORT_TAGS = ( all => \@EXPORT_OK );

=head1 SYNOPSIS

    use Math::Primality::AKS;

    my $n = 123;
    print 'Prime!' if is_aks_prime($n);

=head1 DESCRIPTION

=head1 EXPORT

=head1 FUNCTIONS

=head2 is_aks_prime($n)

Returns 1 if $n is an AKS prime, 0 if it is not.

=cut

sub debug {
    if ( DEBUG or $ENV{DEBUG} ) {
      warn $_[0];
    }
}

sub is_aks_prime($) {
  # http://ece.gmu.edu/courses/ECE746/project/F06_Project_resources/Salembier_Southerington_AKS.pdf
  # http://islab.oregonstate.edu/koc/ece575/04Project2/Halim-Chanleudfa/Report.pdf
  # http://fatphil.org/maths/AKS/
  # http://www.cs.cmu.edu/afs/cs/user/mjs/ftp/thesis-program/2005/rotella.pdf
  # http://www.cse.iitk.ac.in/users/manindra/algebra/primality_v6.pdf

  # This code follows the Rotella 2005 implementation.  Unfortunately that
  # paper was based on the original AKS algorithm, which is *very* slow.  His
  # graphs show hundreds of seconds to test the primality of 3-digit numbers,
  # with a C+GMP implementation.  We therefore expect this Perl+GMPz version
  # to also be completely unusable due to performance, and should be surprised
  # if it proves primality of a 5-digit prime in under an hour.
  #
  # The newer algorithm described in the primality_v6.pdf paper is much faster,
  # though still totally impractical for any actual work without a *lot* of
  # optimization work (see Crandall and Papadopoulos for example).  The
  # variant algorithms of Bernstein become practical in terms of speed, though
  # still do not match the speed of APR-CL or ECPP.

  my $n = GMP->new($_[0]);
  # Step 0 - check that n is a positive integer >= 2;
  if (Rmpz_cmp_ui($n, 2) < 0) {
    debug "fails step 0 of aks - $n is in range\n";
    return 0;  # negative, 0, or 1.
  }
  # Step 1 - check if $n = m^d for some m, d
  # if $n is a power then return 0
  if (Rmpz_perfect_power_p($n)) {
    debug "fails step 1 of aks - $n is a perfect power\n";
    return 0;  # composite
  }

  # This will overestimate r a little, but it's the only good way to do it
  # in GMP without MPFR.  What we want:
  #    floor( log2(n) * log2(n) )
  # What this does:
  #    ceil(log2(n)) * ceil(log2(n))
  # Ex: For 85991 r = 293 instead of 283, with polylimit 290 instead of 275.
  my $logn = Rmpz_sizeinbase($n, 2);
  my $limit = Math::GMPz->new($logn * $logn);

  # Search for first r where order(r, n) > limit
  my $r = Math::GMPz->new(2);
  while (Rmpz_cmp($r, $n) == -1) {
    if(Rmpz_divisible_p($n, $r)) {
        debug "$n is divisible by $r\n";
        return 0;
    }
    my $i = Math::GMPz->new(1);
    my $res = Math::GMPz->new(0);
    my $failed = 0;
    for ( ; Rmpz_cmp($i, $limit) <= 0; Rmpz_add_ui($i, $i, 1)) {
        Rmpz_powm($res, $n, $i, $r);
        if (Rmpz_cmp_ui($res, 1) == 0) {
            $failed = 1;
            last;
        }
    }
    last if !$failed;
    Rmpz_nextprime($r, $r);
  }
  # We've performed trial division for every prime <= r.
  # Hence if n > r*r then n is prime.
  if (Rmpz_cmp($r*$r, $n) >= 0) {
      debug "Found $n is prime while checking for r\n";
      return 1;
  }

    # Polynomial check.  Limit is floor(log2n * sqrt(phi(r))).
    # r is prime hence phi(r) = r-1.
    my $sqrt_phi_r = Math::GMPz->new($r - 1);
    my $plus_one = Rmpz_perfect_square_p($sqrt_phi_r) ? 0 : 1;
    Rmpz_sqrt($sqrt_phi_r, $sqrt_phi_r);
    my $polylimit = $logn * $sqrt_phi_r + $plus_one;

    my $intr = Rmpz_get_ui($r);
    debug "Running poly check, r = $intr  polylimit = $polylimit\n";

    my $final_size = Math::GMPz->new(0);
    Rmpz_mod($final_size, $n, $r);
    my $compare = Math::Primality::BigPolynomial->new();
    $compare->setCoef(Math::GMPz->new(1), $final_size);

    for(my $a = 1; Rmpz_cmp_ui($polylimit, $a) >= 0; $a++) {
        debug "Checking at $a / $polylimit\n";
        my $poly = Math::Primality::BigPolynomial->new( [$a, 1 ]);
        #print "Start:   ", $poly->string, "\n";
        $poly->powmod($n, $n, $intr);
        #print "Finish:  ", $poly->string, "\n";

        $compare->setCoef($a, 0);
        #print "Compare: ", $compare->string, "\n";

        if ( ! $poly->isEqual($compare) ) {
            debug "Found not prime at $a\n";
            return 0;
        }
    }
    return 1;
}

=head1 AUTHORS

Bob Kuo, C<< <bobjkuo at gmail.com> >>
Jonathan "Duke" Leto C<< <jonathan@leto.net> >>

=head1 BUGS

Please report any bugs or feature requests to
C<bug-math-primality-aks at rt.cpan.org>, or through the web interface at
L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Math::Primality::AKS>.
I will be notified, and then you'll automatically be notified of progress on your bug as I make changes.

=head1 THANKS

=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Math::Primality::AKS


You can also look for information at:

=head1 ACKNOWLEDGEMENTS

=head1 COPYRIGHT & LICENSE

Copyright 2009-2010 Bob Kuo, all rights reserved.

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

=cut

exp(0); # End of Math::Primality::AKS
