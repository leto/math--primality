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

Version 0.03_01

=cut

our $VERSION = '0.03_01';

our @EXPORT_OK = qw/is_pseudoprime is_strong_pseudoprime is_strong_lucas_pseudoprime is_prime next_prime prev_prime prime_count/;

our %EXPORT_TAGS = ( all => \@EXPORT_OK );

=head1 SYNOPSIS

    use Math::Primality qw/:all/;

    my $t1 = is_pseudoprime($x,$base);
    my $t2 = is_strong_pseudoprime($x);

    print "Prime!" if is_prime($outrageously_large_prime);

    my $t3 = next_prime($x); 

=head1 DESCRIPTION

Math::Primality implements is_prime() and next_prime() as a replacement for Math::PARI::is_prime().  It uses the GMP library through Math::GMPz.  The is_prime() method is actually a Baillie-PSW primality test which consists of three steps:

=over 4

=item * Check N for small prime divisors p < 1000

=item * Perform a strong Miller-Rabin probable prime test (base 2) on N

=item * Perform a strong Lucas-Selfridge test on N (using the parameters suggested by Selfridge)

=back

At any point the function may return as definitely composite.  If not, N has passed the strong Baillie-PSW test and is either prime or a strong Baillie-PSW pseudoprime.  To date no counterexample (Baillie-PSW strong pseudoprime) is known to exist for N < 10^15.  Baillie-PSW requires O((log n)^3) bit operations.  See L<http://www.trnicely.net/misc/bpsw.html> for a more thorough introduction to the Baillie-PSW test. Also see L<http://mpqs.free.fr/LucasPseudoprimes.pdf> for a more theoretical introduction to the Baillie-PSW test. 

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

=head3 Details

A pseudoprime is a number that satisfies Fermat's Little Theorm, that is, $b^ ($n - 1) = 1 mod $n.

=cut

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

sub debug {
    warn $_[0] if $ENV{DEBUG} or $DEBUG;
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

=head3 Details

A strong pseudoprime to $base is an odd number $n with ($n - 1) = $d * 2^$s that either satisfies

=over 4

=item * $base^$d = 1 mod $n

=item * $base^($d * 2^$r) = -1 mod $n, for $r = 0, 1, ..., $s-1

=back

=head3 Notes

$s and $d are calculated with the helper function _find_s_d() and the second condition is checked by sucessive squaring $base^$d and reducing that mod $n.

=cut

sub is_strong_pseudoprime($;$)
{
    my ($n, $base) = @_;

    $base ||= 2;
    # we should check if we are passed a GMPz object
    $base   = GMP->new($base);
    $n      = GMP->new($n);

    # unnecessary but faster if $n is even
    my $cmp = _check_two_and_even($n);
    return $cmp if $cmp != 2;

    my $m   = _copy($n);
    Rmpz_sub_ui($m,$m,1);

    my ($s,$d) = _find_s_d($m);
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

# given an odd number N find (s, d) such that N = d * 2^s + 1
sub _find_s_d($)
{
    my $m   = GMP->new($_[0]);
    my $s   = Rmpz_scan1($m,1);
    my $d   = GMP->new(0);
    Rmpz_tdiv_q_2exp($d,$m,$s);
    return ($s,$d);
}

=head2 is_strong_lucas_pseudoprime($n)

Returns true if $n is a strong Lucas-Selfridge pseudoprime, false otherwise.  The variable $n should be a Perl
integer or a Math::GMPz object.  

    if ( is_strong_lucas_pseudoprime($n) ) {
        # it's a strong Lucas-Selfridge pseudoprime
    } else {
        # not a strong Lucas-Selfridge psuedoprime
        # i.e. definitely composite
    }

=head3 Details

If we let

=over 4

=item * $D be the first element of the sequence 5, -7, 9, -11, 13, ... for which ($D/$n) = -1.  Let $P = 1 and $Q = (1 - $D) /4

=item * U($P, $Q) and V($P, $Q) be Lucas sequences

=item * $n + 1 = $d * 2^$s + 1

=back

Then a strong Lucas-Selfridge pseudoprime is an odd, non-perfect square number $n with that satisfies either

=over 4

=item * U_$d = 0 mod $n

=item * V_($d * 2^$r) = 0 mod $n, for $r = 0, 1, ..., $s-1

=back

=head3 Notes

($d/$n) refers to the Legendre symbol.
The tuple ($D, $P, $Q) is determined by the helper function _find_dpq_selfridge(). 
$d and $s are determined by the helper function _find_s_d().

=cut

sub is_strong_lucas_pseudoprime($)
{
    my ($n) = @_;
    $n      = GMP->new($n);
    # we also need to weed out all N < 3 and all even N 
    my $cmp = _check_two_and_even($n);
    return $cmp if $cmp != 2;
    # weed out all perfect squares
    if ( Rmpz_perfect_square_p($n) ) {
        return 0;
    }
    # determine Selfridge parameters D, P and Q
    my ($D, $P, $Q) = _find_dpq_selfridge($n);
    if ($D == 0) {  #_find_dpq_selfridge found a factor of N
      return 0;
    }
    my $m = _copy($n);
    Rmpz_add_ui($m, $m, 1);

    # determine s and d such that m = d * 2^s + 1
    my ($s,$d) = _find_s_d($m);
    # compute U_d and V_d
    # initalize U, V, U_2m, V_2m
    my $U = GMP->new(1);     # U = U_1 = 1
    my $V = GMP->new($P);    # V = V_1 = P
    my $U_2m = GMP->new(1);  # U_2m = U_1
    my $V_2m = GMP->new($P); # V_2m = P
    # initalize Q values (eventually need to calculate Q^d, which will be used in later stages of test)
    my $Q_m = GMP->new($Q);
    my $Q2_m = GMP->new(2 * $Q);  # Really 2Q_m, but perl will barf with a variable named like that
    my $Qkd = GMP->new($Q);
    # start doubling the indicies!
    my $dbits = Rmpz_sizeinbase($d,2);
    for (my $i = 1; $i < $dbits; $i++) {  #since d is odd, the zeroth bit is on so we skip it
      # U_2m = U_m * V_m (mod N)
      Rmpz_mul($U_2m, $U_2m, $V_2m);  # U_2m = U_m * V_m
      Rmpz_mod($U_2m, $U_2m, $n);     # U_2m = U_2m mod N
      # V_2m = V_m * V_m - 2 * Q^m (mod N)
      Rmpz_mul($V_2m, $V_2m, $V_2m);  # V_2m = V_2m * V_2m
      Rmpz_sub($V_2m, $V_2m, $Q2_m);  # V_2m = V_2m - 2Q_m
      Rmpz_mod($V_2m, $V_2m, $n);     # V_2m = V_2m mod N
      # calculate powers of Q for V_2m and Q^d (used later)
      # 2Q_m = 2 * Q_m * Q_m (mod N)
      Rmpz_mul($Q_m, $Q_m, $Q_m);     # Q_m = Q_m * Q_m
      Rmpz_mod($Q_m, $Q_m, $n);       # Q_m = Q_m mod N
      Rmpz_mul_2exp($Q2_m, $Q_m, 1);  # 2Q_m = Q_m * 2
      if (Rmpz_tstbit($d, $i)) {       # if bit i of d is set
        # add some indicies
        # initalize some temporary variables
        my $T1 = GMP->new(0);
        my $T2 = GMP->new(0);
        my $T3 = GMP->new(0);
        my $T4 = GMP->new(0);
        # this is how we do it
        # U_(m+n) = (U_m * V_n + U_n * V_m) / 2
        # V_(m+n) = (V_m * V_n + D * U_m * U_n) / 2
        Rmpz_mul($T1, $U_2m, $V);     # T1 = U_2m * V
        Rmpz_mul($T2, $U, $V_2m);     # T2 = U * V_2m
        Rmpz_mul($T3, $V_2m, $V);     # T3 = V_2m * V
        Rmpz_mul($T4, $U_2m, $U);     # T4 = U_2m * U
        Rmpz_mul_si($T4, $T4, $D);    # T4 = T4 * D = U_2m * U * D
        Rmpz_add($U, $T1, $T2);       # U = T1 + T2 = U_2m * V - U * V_2m
        if (Rmpz_odd_p($U)) {         # if U is odd
          Rmpz_add($U, $U, $n);       # U = U + n
        }
        Rmpz_fdiv_q_2exp($U, $U, 1);  # U = floor(U / 2)
        Rmpz_add($V, $T3, $T4);       # V = T3 + T4 = V_2m * V + U_2m * U * D
        if (Rmpz_odd_p($V)) {         # if V is odd
          Rmpz_add($V, $V, $n);       # V = V + n 
        }
        Rmpz_fdiv_q_2exp($V, $V, 1);  # V = floor(V / 2)
        Rmpz_mod($U, $U, $n);         # U = U mod N
        Rmpz_mod($V, $V, $n);         # V = V mod N
        # Get our Q^d calculating on (to be used later)
        Rmpz_mul($Qkd, $Qkd, $Q_m);   # Qkd = Qkd * Q_m
        Rmpz_mod($Qkd, $Qkd, $n);     # Qkd = Qkd mod N
      }
    }
    # if U_d or V_d = 0 mod N, then N is prime or a strong Lucas pseudoprime
    if(Rmpz_sgn($U) == 0 || Rmpz_sgn($V) == 0) {
      return 1;
    }
    # ok, if we're still here, we have to compute V_2d, V_4d, V_8d, ..., V_{2^(s-1)*d}
    # initalize 2Qkd
    my $Q2kd = GMP->new;
    Rmpz_mul_2exp($Q2kd, $Qkd, 1);    # 2Qkd = 2 * Qkd
    # V_2m = V_m * V_m - 2 * Q^m (mod N)
    for (my $r = 1; $r < $s; $r++) {
      Rmpz_mul($V, $V, $V);      # V = V * V;
      Rmpz_sub($V, $V, $Q2kd);   # V = V - 2Qkd
      Rmpz_mod($V, $V, $n);      # V = V mod N
      # if V = 0 mod N then N is a prime or a strong Lucas pseudoprime
      if(Rmpz_sgn($V) == 0) {
        return 1;
      }
      # calculate Q ^(d * 2^r) for next r (unless on final iteration)
      if ($r < ($s - 1)) {
        Rmpz_mul($Qkd, $Qkd, $Qkd);     # Qkd = Qkd * Qkd
        Rmpz_mod($Qkd, $Qkd, $n);       # Qkd = Qkd mod N
        Rmpz_mul_2exp($Q2kd, $Qkd, 1);  # 2Qkd = 2 * Qkd
      } 
    }
    # otherwise N is definitely composite 
    return 0;
}

# selfridge's method for finding the tuple (D,P,Q) for is_strong_lucas_pseudoprime
sub _find_dpq_selfridge($) {
  my $n = GMP->new($_[0]);
  my ($d,$sign,$wd) = (5,1,0);
  my $gcd = Math::GMPz->new;

  # determine D
  while (1) {
    $wd = $d * $sign;

    Rmpz_gcd_ui($gcd, $n, abs $wd);
    if ($gcd > 1 && Rmpz_cmp($n, $gcd) > 0) {
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
    # need code to make sure we don't overflow $d, may never actually happen
    ### TODO ###
  }
  # P = 1
  my ($p,$q) = (1,0);
  {
      use integer;
      # Q = (1 - D) / 4
      $q = (1 - $wd) / 4;
  }
  debug "found P and Q: ($p, $q)";
  return ($wd, $p, $q);
}

# method returns 0 if N < two or even, returns 1 if N == 2
# returns 2 if N > 2 and odd
sub _check_two_and_even($) {
  my $n = GMP->new($_[0]);

  my $cmp = Rmpz_cmp_ui($n, 2 );
  return 1 if $cmp == 0;
  return 0 if $cmp < 0;
  return 0 if Rmpz_even_p($n);
  return 2;
}

=head2 is_prime($n)

Returns true if number is prime, false if number is composite.

    if ( is_prime($n) ) {
        # it's a prime
    } else {
        # definitely composite
    }

=head3 Details

is_prime() is implemented using the BPSW algorithim which is a combination of two probable-prime 
algorithims, the strong Miller-Rabin test and the strong Lucas-Selfridge test.  While no
psuedoprime has been found for N < 10^15, this does not mean there is not a pseudoprime.

=head3 Notes

The strong Miller-Rabin test is implemented by is_strong_pseudoprime(). The strong Lucas-Selfridge test is implemented
by is_strong_lucas_pseudoprime().

=cut

sub is_prime($) {
  my $n = GMP->new($_[0]);
  # TODO: 
  # trial division of n up to some small number (perhaps a thousand)

  # the lucas test is stronger so do it first
  return is_strong_lucas_pseudoprime($n) && is_strong_pseudoprime($n,2);
}

=head2 next_prime($x)

Given a number, produces the next prime number.

    my $q = next_prime($n);

=head3 Details

Each next greatest odd number is checked until one is found to be prime

=head3 Notes

Checking of primality is implemented by is_prime()

=cut

sub next_prime($) {
  my $n = GMP->new($_[0]);
  my $cmp = Rmpz_cmp_ui($n, 2 ); #check if $n < 2
  if ($cmp < 0) {
    return GMP->new(2);
  }
  if (Rmpz_odd_p($n)) {         # if N is odd
    Rmpz_add_ui($n, $n, 2);     # N = N + 2
  } else {
    Rmpz_add_ui($n, $n, 1);     # N = N + 1
  }
  # N is now the next odd number
  while (1) {
    return $n if is_prime($n);  # check primality of that number, return if prime
    Rmpz_add_ui($n, $n, 2);     # N = N + 2
  }
}

=head2 prev_prime($x)

Given a number, produces the previous prime number.

    my $q = prev_prime($n);

=head3 Details

Each previous odd number is checked until one is found to be prime.  prev_prime(2) or for any number less than 2 returns undef

=head3 Notes

Checking of primality is implemented by is_prime()

=cut

sub prev_prime($) {
  my $n = GMP->new($_[0]);
  my $cmp = Rmpz_cmp_ui($n, 3);   # compare N with 3
  if ($cmp == 0) {                # N = 0
    return GMP->new(2);
  } elsif ($cmp < 0) {            # N < 0
    return undef;
  } else {
    if (Rmpz_odd_p($n)) {         # if N is odd
      Rmpz_sub_ui($n, $n, 2);     # N = N - 2
    } else {
      Rmpz_sub_ui($n, $n, 1);     # N = N - 1
    }
    # N is now the previous odd number
    while (1) {
      return $n if is_prime($n);  # check primality of that number, return if prime
      Rmpz_sub_ui($n, $n, 2);     # N = N - 2
    } 
  }
}

=head2 prime_count($n)

Returns the count of the number of primes less than or equal to $n. This is the
prime counting function.

=cut

sub prime_count($) {
  my $n = GMP->new($_[0]); # check if $n needs to be a Math::GMPz object
  my $primes = 0;
  return 0 if $n <= 1;

  for (my $i = GMP->new(0); Rmpz_cmp($i, $n) <= 0; Rmpz_add_ui($i, $i, 1)) {
    $primes++ if is_prime($i);
  }
  return $primes;
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
