use Test::More tests => 9;

use Math::Brent qw(Brentzero);
use Math::Utils qw(:compare :polynomial);
use strict;
use warnings;

my($x, $y);
my $fltcmp = generate_fltcmp(5e-7);
my $brent_tol = 1e-8;

sub wobble
{
	my($t) = @_;
	return $t - cos($t);
}

($x, $y) = Brentzero(0, 1, \&wobble, $brent_tol);
ok(&$fltcmp($x, 0.739085133) == 0, "Wobble(), ($x, $y)");

#
# Some simple polynomials.
#
my $eqn1 = sub {my($x) = @_; return pl_evaluate(-12, -11, 2, 1);};
my $eqn2 = sub {my($x) = @_; return pl_evaluate(-1, -2, 11, 12);};
my $eqn3 = sub {my($x) = @_; return pl_evaluate(11, 10, 8, 5, 1);};

#
# Roots of the cubic eqn1.
#
($x, $y) = Brentzero(1, 5, $eqn1);
ok(&$fltcmp($x, 3.0) == 0, "Anon sub 1, ($x, $y)");

($x, $y) = Brentzero(-2, 0, $eqn1);
ok(&$fltcmp($x, -1.0) == 0, "Anon sub 1, ($x, $y)");

($x, $y) = Brentzero(-5, -3.5, $eqn1);
ok(&$fltcmp($x, -4.0) == 0, "Anon sub 1, ($x, $y)");

#
# Roots of the cubic eqn2.
#
($x, $y) = Brentzero(0.0, 0.5, $eqn2);
ok(&$fltcmp($x, 0.33333333) == 0, "Anon sub 2, ($x, $y)");

($x, $y) = Brentzero(-0.75, 0.25, $eqn2);
ok(&$fltcmp($x, -0.25) == 0, "Anon sub 2, ($x, $y)");

($x, $y) = Brentzero(-3, -0.5, $eqn2);
ok(&$fltcmp($x, -1) == 0, "Anon sub 2, ($x, $y)");

#
# First root of eqn3...
#
($x, $y) = Brentzero(-5, -2, $eqn3);
ok(&$fltcmp($x, -3.08054627) == 0, "Anon sub 3, ($x, $y)");

#
# ... and now its second root (the other two are complex).
#
($x, $y) = Brentzero(-2, 0, $eqn3);
ok(&$fltcmp($x, -3.08054627) == 0, "Anon sub 3, ($x, $y)");

1;
