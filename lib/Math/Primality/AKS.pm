package Math::Primality::AKS;
use warnings;
use strict;

use Math::GMPz qw/:mpz/;
use base 'Exporter';
use Carp qw/croak/;

#use Math::Primality;

use constant DEBUG => 0;

use constant GMP => 'Math::GMPz';

=head1 NAME

Math::Primality::AKS - The Agrawal-Kayal-Saxena primality testing algorithm

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';
$VERSION = eval $VERSION;

our @EXPORT_OK = qw/is_aks_prime/;

our %EXPORT_TAGS = ( all => \@EXPORT_OK );

=head1 SYNOPSIS

    use Math::Primality::AKS;

    my $n = 123;
    print 'Prime!' if is_aks_prime($n);

=head1 DESCRIPTION

=head1 EXPORT

=head1 FUNCTIONS

=head2 aks($n)

=cut

sub debug {
    if ( DEBUG ) {
      warn $_[0];
    }
}

sub is_aks_prime($) {
  # http://ece.gmu.edu/courses/ECE746/project/F06_Project_resources/Salembier_Southerington_AKS.pdf
  # http://islab.oregonstate.edu/koc/ece575/04Project2/Halim-Chanleudfa/Report.pdf
  # http://fatphil.org/maths/AKS/
  my $n = GMP->new($_[0]);
  # Step 1 - check if $n = m^d for some m, d
  # if $n is a power then return 0
  if (Rmpz_perfect_power_p($n)) {
    debug "fails step 1 of aks - $n is a perfect power" if DEBUG;
    return 0;  # composite
  }
  # Step 2 - find smallest $r > (log ($n))^2 such that (x + a) ^ $n = x^$n + a (mod x^$r - 1, n)
  my $r = GMP->new(_find_smallest_r($n));
  # Step 3 - if gcd($a, $n) != 1 for all 1<= $a <= $r then return 0
  for (my $a = GMP->new(1); Rmpz_cmp($a, $r) <= 0; Rmpz_add_ui($a, $a, 1)) {
    my $gcd_rslt = GMP->new();
    Rmpz_gcd($gcd_rslt, $a, $n);  # $gcd_rslt = gcd($a, $n)
    if (Rmpz_cmp_ui($gcd_rslt, 1) > 0 && Rmpz_cmp($gcd_rslt, $n) < 0) {  # if $gcd_rslt > 1 && $gcd_rslt < $n
      debug "fails step 3 of aks - $gcd_rslt is a factor of $n" if DEBUG;
      return 0; # composite
    } 
  }
  # Step 4 - if $n <= $r then return 0
  if (Rmpz_cmp($n, $r) <= 0) {
    debug "step 4 of aks - $n is prime" if DEBUG;
    return 2; # prime
  }
  # Step 5 - for $a = 1 to floor (2*sqrt(totient($r))*log($n)) do
    # if ( (x + a) ^ $n = x^$n + a (mod x^$r - 1, $n)) then return 0
  # compute $max_a = 2 * sqrt(totient($r)) * log($n)
  my $max_a = GMP->new(_calculate_max_a($r, $n));
  my $a = GMP->new(1);
  for (; Rmpz_cmp($a, $max_a) <= 0; Rmpz_add_ui($a, $a, 1)) {
    # if ( (x + a) ^ $n = x^$n + a (mod x^$r - 1, $n)) then return 0
    if (!_poly_eq_holds($a, $n, $r)) {
      debug "Step 5 detected composite" if DEBUG;
      return 0;
    }
  } 
  # Step 6 - return 2
  return 2; # prime
}

sub _calculate_max_a($;$) {
  my $r = $_[0];
  my $n = $_[1];
  # since $r is prime, totient($r) is just $r - 1
  my $totient_r = GMP->new();
  Rmpz_sub_ui($totient_r, $r, 1);  # $totient_r = $r - 1
  # compute ceil(log ($n))
  my $logn = GMP->new(_Rmpz_logbase2cl($n)); # $logn = ceil(log ($n))
  # compute 2 * ceil(log ($n))
  my $logn2 = GMP->new(); 
  Rmpz_mul_ui($logn2, $logn, 2);   # logn2 = $logn * 2
  # to save us from using floating point, calculate totient_($r) * (2 * log ($n))^2 and then take the sqrt
  Rmpz_mul($logn2, $logn2, $logn2); 
  # compute 4 * totient($r) * (log ($n))^2
  my $squared = GMP->new();
  Rmpz_mul($squared, $totient_r, $logn2); 
  # compute 2 * sqrt(totient($r)) * log ($n)
  my $max_a = GMP->new();
  Rmpz_sqrt($max_a, $squared);
  debug "max_a is $max_a" if DEBUG;
  return $max_a;
}

sub _poly_eq_holds($) {
  my $a = $_[0];
  my $n = $_[1];
  my $r = $_[2];
  # my @r_mod = _sparse_polynomial($r, -1);
  # my @x_a = [$a, 1];
  # my @rhs = pow_mod(\@x_a, $n, \@r_mod);
  # my @x_n = _sparse_polynomial($n, $a);
  # my @lhs = mod(\@x_n, \@r_mod);
  # if (nmod(\@lhs, $n) == nmod(\$rhs, $n)) {
  # do something significant
  # }
  return 2;
}

sub _find_smallest_r($) {
  # try out successive values of $r and test if $n^$k != 1 mod $r for every $k <= 4(log n)^2
  my $n = $_[0];
  # find $max_k
  my $logn = GMP->new(_Rmpz_logbase2cl($n)); # $logn = log $n (base 2)
  my $logn_sqrd = GMP->new();
  Rmpz_mul($logn_sqrd, $logn, $logn);        # $logn_sqrd = $logn * $logn
  my $max_k = GMP->new();
  Rmpz_mul_si($max_k, $logn_sqrd, 4);        # $max_k = 4 * $logn_sqrd
  debug "max_k is $max_k" if DEBUG;
  my $r = GMP->new(2);
  while (1) {
    if (_r_good($r, $max_k, $n)) {
      last;
    }
    Rmpz_add_ui($r, $r, 1);                  # $r = $r + 1
  }
  debug "r is $r" if DEBUG;
  return $r;
}

# private function _r_good ($r, $max_k, $n)
# return 0 is this $r isn't good for a certain $n, 1 if it is
sub _r_good($;$;$) {
  my $r = $_[0];
  my $max_k = $_[1];
  my $n = $_[2];

  my $k = GMP->new(1);
  my $mod_rslt = GMP->new();

  while (Rmpz_cmp($k, $max_k) <= 0) {
    Rmpz_powm($mod_rslt, $n, $k, $r);       # $mod_rslt = $n ^ $k mod $r
    if (Rmpz_cmp_ui($mod_rslt, 1) == 0) {   # if $mod_rslt == 0
      return 0;
    }
    Rmpz_add_ui($k, $k, 1);                 # $k = $k + 1
  }
  return 1;
}

# thanks GMP for not having a log function...
# this nonsense with the subtracting is for the ceiling (cl) and floor (fl)
sub _Rmpz_logbase2cl($) {
  my $n = $_[0];
  my ($double, $si) = Rmpz_get_d_2exp($n);  # $double * 2^$si ~= $op (with 0.5 <= abs($double) < 1)
  if ($double == 0.5) {
    $si--;
  } 
  return $si;
}

sub _Rmpz_logbase2fl($) {
  my $n = $_[0];
  my ($double, $si) = Rmpz_get_d_2exp($n);  # $double * 2^$si ~= $op (with 0.5 <= abs($double) < 1)
  if ($double == 0.5) {
    $si++;
  } else {
    $si--;
  }
  return $si;
}

=head1 AUTHORS

Bob Kuo, C<< <bobjkuo at gmail.com> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-math-primality-aks at rt.cpan.org>,
or through the web interface at
L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Math::Primality::AKS>.  I will be 
notified, and then you'll automatically be notified of progress on your bug as I
make changes.

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
