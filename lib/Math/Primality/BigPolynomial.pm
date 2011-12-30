package Math::Primality::BigPolynomial;

use strict;
use warnings;
use Math::GMPz qw/:mpz/;

sub new {
    my $self              = {};
    my $class             = shift;
    my $construction_junk = shift;
    if ($construction_junk) {
        if ( ref($construction_junk) eq 'ARRAY' ) {
            $self->{COEF}   = $construction_junk;
            $self->{DEGREE} = scalar(@$construction_junk);
        } else {
            $self->{DEGREE} = $construction_junk;
            my @a = [];
            for ( my $i = 0 ; $i < $construction_junk ; $i++ ) {
                push @a, Math::GMPz->new(0);
            }
            $self->{COEF} = \@a;
        }
    }
    else {
        $self->{COEF}   = [ Math::GMPz->new(0) ];
        $self->{DEGREE} = 1;
    }
    bless( $self, $class );
    return $self;
}

sub coef {
    my $self = shift;
    if (@_) { @{ $self->{COEF} } = @_ }
    return @{ $self->{COEF} };
}

sub degree {
    my $self = shift;
    if (@_) { $self->{DEGREE} = shift }
    return $self->{DEGREE};
}

sub getCoef {
    my $self = shift;
    my $i    = shift;
    if ( $i > $self->degree() ) {
        return 0;
    }
    return ${ $self->{COEF} }[$i];
}

sub isEqual {
    my $self             = shift;
    my $other_polynomial = shift;
    if ( $self->degree() != $other_polynomial->degree() ) {
        return 0;
    }
    for ( my $i = 0 ; $i < $self->degree() ; $i++ ) {
        if ( $self->getCoef($i) != $other_polynomial->getCoef($i) ) {
            return 0;
        }
    }
    return 1;
}

sub setCoef {
    my $self     = shift;
    my $new_coef = shift;
    my $index    = shift;
    if ( $index < 0 ) {
        die "coef is less than 0";
    }

    if ( $index > $self->degree() ) {
        for ( my $j = $self->degree() + 1 ; $j <= $index ; $j++ ) {
            push @{ $self->{COEF} }, Math::GMPz->new(0);
        }
        push @{ $self->{COEF} }, $new_coef;
        $self->degree($index);
    }
    else {
        ${ $self->{COEF} }[$index] = $new_coef;
    }
}

sub compact {
    my $self = shift;
    my $i    = 0;
  LOOP: for ( $i = $self->degree() - 1 ; $i > 0 ; $i-- ) {
        if ( Math::GMPz::Rmpz_cmp_ui( $self->getCoef($i), 0 ) != 0 ) {
            last LOOP;
        }
        pop @{ $self->{COEF} };
    }
    if ( $i != $self->degree() ) {
        $self->degree( $i + 1 );
    }
}

sub clear {
    my $self = shift;
    $self->{COEF}   = undef;
    $self->{DEGREE} = undef;
    $self->{COEF}   = [ Math::GMPz->new(0) ];
    $self->{DEGREE} = 1;
}

sub mpz_poly_mod_mult {
    my $self = shift;
    my ( $rop, $x, $y, $mod, $polymod ) = @_;

    $rop->clear();

    my $xdeg   = $x ? $x->degree() : 0;
    my $ydeg   = $y ? $y->degree() : 0;
    my $maxdeg = $xdeg < $ydeg ? $ydeg : $xdeg;

  LOOP: for ( my $i = 0 ; $i < $polymod ; $i++ ) {
        my $sum  = Math::GMPz->new(0);
        my $temp = Math::GMPz->new(0);
        for ( my $j = 0 ; $j <= $i ; $j++ ) {
            Rmpz_add( $temp,
                $y->getCoef( $i - $j ) + $y->getCoef( $i + $polymod - $j ) );
            Rmpz_mul( $temp, $x->getCoef($j), $temp );
            Rmpz_add( $sum, $sum, $temp );
        }

        for ( my $j = 0 ; $j < ( $i + $polymod ) ; $j++ ) {
            Rmpz_mul( $temp, $x->getCoef($j),
                $y->getCoef( $i + $polymod - $j ) );
            Rmpz_add( $sum, $sum, $temp );
        }

        Rmpz_mod( $temp, $sum, $mod );
        $rop->setCoef( $temp, $i );

        if ( $i > $maxdeg && Rmpz_cmp_ui( $sum, 0 ) == 0 ) {
            last LOOP;
        }
    }

    $rop->compact();
}

sub mpz_poly_mod_power {
    my $self = shift;
    my ( $rop, $x, $power, $mult_mod, $poly_mod ) = @_;

    $rop->clear();
    $rop->setCoef( Math::GMPz->new(1), 0 );

    my $i = Rmpz_sizeinbase( $power, 2 );

  LOOP: for ( ; $i >= 0 ; $i-- ) {
        mpz_poly_mod_mult( $rop, $rop, $rop, $mult_mod, $poly_mod );

        if ( Rmpz_tstbit( $power, $i ) ) {
            mpz_poly_mod_mul( $rop, $rop, $x, $mult_mod, $poly_mod );
        }

        if ( $i == 0 ) {
            last LOOP;
        }
    }

    $rop->compact();
}

1;
