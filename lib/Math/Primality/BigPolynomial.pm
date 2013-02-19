package Math::Primality::BigPolynomial;

use strict;
use warnings;
use Math::GMPz qw/:mpz/;

# ABSTRACT: Big Polynomials

=head1 NAME

Math::Primality::BigPolynomials - Polynomials with BigInts

=cut


sub new {
    my $self   = {};
    my $class  = shift;
    my $init   = shift;

    $self->{ZERO} = Math::GMPz->new(0);

    if ($init) {
        my $type = ref $init;
        if ( $type eq 'ARRAY' ) {
            die "Initialization array must be non-empty" unless @$init > 1;
            $self->{COEF} = $init;
        } elsif ($type eq 'Math::Primality::BigPolynomial') {
            $self->{COEF} = [ map { Math::GMPz->new($_) } $init->coef ];
        } else {
            my $a = [];
            for ( my $i = 0 ; $i < $init ; $i++ ) {
                push @$a, $self->{ZERO};
            }
            $self->{COEF} = $a;
        }
    }
    else {
        $self->{COEF} = [ $self->{ZERO} ];
    }
    bless( $self, $class );
    return $self;
}

sub copy {
  my $self = shift;
  return Math::Primality::BigPolynomial->new($self);
}

sub string {
  my $self = shift;
  my @coefs = $self->coef;
  my $string = '';
  foreach my $i (reverse 1 .. $#coefs) {
    my $c = $coefs[$i];
    next unless $c;
    $string .= $c if $c != 1;
    $string .= ($i > 1) ? "x^$i" : "x";
    $string .= ' + ';
  }
  $string .= $coefs[0];
  return $string;
}

sub coef {
    my $self = shift;
    if (@_) { @{ $self->{COEF} } = @_ }
    return @{ $self->{COEF} };
}

sub degree {
    my $self = shift;
    return (scalar @{$self->{COEF}} - 1);
}

sub getCoef {
    my $self = shift;
    my $i    = shift;
    return if !defined $i || $i < 0;
    return $self->{ZERO} if $i >= scalar @{$self->{COEF}};
    return $self->{COEF}->[$i];
}

sub isEqual {
    my $self  = shift;
    my $other = shift;
    return 0 unless $self->degree == $other->degree;
    foreach my $i (0 .. $self->degree) {
      return 0 unless $self->getCoef($i) == $other->getCoef($i);
    }
    return 1;
}

sub setCoef {
    my $self     = shift;
    my $new_coef = shift;
    my $index    = shift;
    die "setCoef: coef not given" unless defined $index;
    die "setCoef: coef $index is negative" if $index < 0;

    for ( my $j = $self->degree() + 1 ; $j < $index ; $j++ ) {
        $self->{COEF}->[$j] = $self->{ZERO};
    }
    $self->{COEF}->[$index] = Math::GMPz->new($new_coef);
}

sub compact {
    my $self = shift;
    while (scalar @{$self->{COEF}} > 1 && $self->{COEF}->[-1] == 0) {
      pop @{ $self->{COEF} };
    }
    return $self;
}

sub clear {
    my $self = shift;
    $self->{COEF}   = [ $self->{ZERO} ];
}

sub mulmod {
    my ($self, $copy_y, $mod, $polymod ) = @_;
    die "mulmod: first argument must be a poly!" unless ref($copy_y) eq 'Math::Primality::BigPolynomial';
    die "mulmod: mod must be defined and > 0!" unless $mod;
    die "mulmod: polymod must be defined and > 0!" unless $polymod;

    # Bypassing getCoef is a 2-3x speedup.
    my @x = @{$self->{COEF}};
    my @y = @{$copy_y->{COEF}};
    my $xdeg   = $#x;
    my $ydeg   = $#y;
    my $maxdeg = $xdeg < $ydeg ? $ydeg : $xdeg;
    # Track which coefficients have non-zero components.
    my @xset = map { Rmpz_sgn($_) > 0 } @x;
    my @yset = map { Rmpz_sgn($_) > 0 } @y;

    my @res = map { Math::GMPz->new(0) } 0 .. $polymod-1;

    for (my $ix = 0; $ix <= $xdeg; $ix++) {
      next unless $xset[$ix];
      for (my $iy = 0; $iy <= $ydeg; $iy++) {
        next unless $yset[$iy];
        Rmpz_addmul( $res[ ($ix + $iy) % $polymod ], $x[$ix], $y[$iy] );
      }
    }

    $self->{COEF} = [ map { $_ % $mod } @res ];
    $self->compact();
}

sub powmod {
  my ($self, $power, $mult_mod, $poly_mod) = @_;
  die "mpz_poly_mod_power: polymod must be defined!" unless $poly_mod;

  my $x = Math::Primality::BigPolynomial->new($self);
  $self->clear();
  $self->setCoef(1, 0);

  my $p = Math::GMPz->new($power);
  while (Rmpz_sgn($p) > 0) {
    $self->mulmod($x, $mult_mod, $poly_mod) if Rmpz_odd_p($p);
    Rmpz_div_2exp($p, $p, 1);
    $x->mulmod($x, $mult_mod, $poly_mod) if Rmpz_sgn($p) > 0;
  }
  $self->compact();
}

1;
