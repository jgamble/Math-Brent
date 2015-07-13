# $Id: Brent.pm,v 1.1 1995/12/26 10:06:36 willijar Exp $ 
=head1 NAME

  Math::Brent - Single Dimensional Function Minimisation

=head1 SYNOPSIS

    use Math::Brent qw(FindMinima BracketMinimum Brent Minimise1D);

    my ($x, $y) = Minimise1D($guess, $scale, \&func, $tol, $itmax);
    my ($ax, $bx, $cx, $fa, $fb, $fc) = BracketMinimum($ax, $bx, \&func);
    my ($x, $y) = Brent($ax, $bx, $cx, \&func, $tol, $itmax);

=head1 DESCRIPTION

This is an implementation of Brent's method for One-Dimensional
minimisation of a function without using derivatives. This algorithm
cleverly uses both the Golden Section Search and parabolic
interpolation.

=head2 FUNCTIONS

The functions may be imported by name, or by using the export
tag "all".

=head3 Brent()

Given a function reference B<\&func> and a
bracketing triplet of abcissas B<$ax>, B<$bx>, B<$cx> (such that
B<$bx> is between B<$ax> and B<$cx> and B<func($bx)> is less than both
B<func($ax)> and B<func($cx)>), Brent() isolates the minimum to a fractional
precision of about B<$tol> using Brent's method. A maximum number of
iterations B<$itmax> may be specified for this search - it defaults to
100. Returned is an array consisting of the abcissa of the minum and
the function value there.

=head3 BracketMinimum()

Given a function reference B<\&func> and
distinct initial points B<$ax> and B<$bx> searches in the downhill
direction (defined by the function as evaluated at the initial points)
and returns an array of the three points B<$ax>, B<$bx>, B<$cx> which
bracket the minimum of the function and the function values at those
points.

=head3 Minimise1D()

Provides a simple interface to the above
two routines. Given a function B<\&func>, an initial guess for its
minimum, and its scaling (B<$guess>, B<$scale>) this routine isolates
the minimum to a fractional precision of about B<$tol> using Brent's
method. A maximum number of iterations B<$itmax> may be specified for
this search - it defaults to 100. It returns an array consisting of
the abcissa of the minum and the function value there.

=head1 EXAMPLE

    use Math::Brent qw(Minimise1D);

    sub sinc {
      my $x = shift ;
      return $x ? sin($x)/$x: 1;
    }

    my($x, $y) = Minimise1D(1, 1, \&sinc, 1e-7);
    print "Minimum of sinc($x) = $y\n";

produces the output

    Minimum is func(5.236068) = -.165388470697432

Anonymous subroutines may also be used as the function reference:

    my $cubic_ref = sub {my($x) = @_; return 6.25 + $x*$x*(-24 + $x*8));};

    my($x, $y) = Minimise1D(3, 1, $cubic_ref);
    print "Minimum of the cubic at $x = $y\n";
    
=head1 BUGS

Let me know of any problems.

=head1 AUTHOR

John A.R. Williams B<J.A.R.Williams@aston.ac.uk>

John M. Gamble B<jgamble@cpan.org> (current maintainer)

=head1 SEE ALSO

"Numerical Recipies: The Art of Scientific Computing"
W.H. Press, B.P. Flannery, S.A. Teukolsky, W.T. Vetterling.
Cambridge University Press. ISBN 0 521 30811 9.

=cut

package Math::Brent;

use strict;
use warnings;
use 5.8.3;

use Exporter;
our (@ISA, @EXPORT_OK, %EXPORT_TAGS);
@ISA = qw(Exporter);
%EXPORT_TAGS = (
	all => [qw(
		FindMinima
		BracketMinimum
		Brent Minimise1D
	) ],
);

@EXPORT_OK = ( @{ $EXPORT_TAGS{all} } );

our $VERSION = 0.04;

use Math::VecStat qw(max min);
use Math::Fortran qw(sign);
use Carp;

sub Minimise1D
{
    my ($guess, $scale, $func, $tol, $itmax) = @_;
    my ($a, $b, $c) = BracketMinimum($guess - $scale, $guess + $scale, $func);

    return Brent($a, $b, $c, $func, $tol, $itmax);
}

#----------------------------------------------------------------------
# BracketMinimum
#
# BracketMinimum is MNBRAK minimum bracketing routine from section 10.1
# of Numerical Recipies
#
# Given a function func, and distinct initial points ax & bx this
# routine searches in the downhill direction and returns new points ax,
# bx, cx which bracket the minimum. The function values at the 3 points
# are returned in fa, fb, fc respectively.
#
my $GOLD= 1.618034; # default magnification ratio for intervals
my $GLIMIT= 100.0; # Max magnification for parabolic fit step
my $TINY= 1E-20;

sub BracketMinimum
{
    my ($ax, $bx, $func) = @_;
    my ($fa, $fb) = (&$func($ax), &$func($bx));

    #
    # Swap the a and b values if we weren't going in
    # a downhill direction.
    #
    if ($fb > $fa)
    {
	my $t = $ax; $ax = $bx; $bx = $t; $t = $fa; $fa = $fb; $fb = $t;
    }

    my $cx = $bx + $GOLD * ($bx - $ax);
    my $fc = &$func($cx);

	# Loop here until we bracket
    while ($fb >= $fc)
    {
	#
	# Compute U by parabolic extrapolation from
	# a, b, c. TINY used to prevent div by zero
	#
	my $r = ($bx - $ax) * ($fb - $fc);
	my $q = ($bx - $cx) * ($fb - $fa);
	my $u = $bx - (($bx - $cx) * $q - ($bx - $ax) * $r)/
	    (2.0 * sign(max(abs($q - $r), $TINY), $q - $r));

	my $ulim = $bx + $GLIMIT * ($cx - $bx); # We won't go further than this
	my $fu;

	# parabolic U between B & C - try it
	if (($bx-$u) * ($u-$cx) > 0.0)
	{
	    $fu = &$func($u);

	    if ($fu < $fc)
	    {
		# Minimum between B & C
		$ax = $bx; $fa = $fb; $bx = $u;  $fb = $fu;
		next;
	    }
	    elsif ($fu > $fb)
	    {
		# Minimum between A & U
		$cx = $u; $fc = $fu;
		next;
	    }

	    $u = $cx + $GOLD * ($cx - $bx);
	    $fu = &$func($u);
	}
	elsif (($cx - $u) * ($u - $ulim) > 0)
	{
		# parabolic  fit between C and limit
	    $fu = &$func($u);

	    if ($fu < $fc)
	    {
		$bx = $cx; $cx = $u;
		$u = $cx + $GOLD*($cx-$bx);
		$fb = $fc; $fc = $fu;
		$fu = &$func($u);
	    }
	}
	elsif (($u - $ulim) * ($ulim - $cx) >= 0)
	{
	    # Limit parabolic U to maximum
	    $u = $ulim;
	    $fu = &$func($u);
	}
	else
	{
	    # eject  parabolic U, use default magnification
	    $u = $cx + $GOLD * ($cx - $bx);
	    $fu = &$func($u);
	}

	# Eliminate oldest point & continue
	$ax = $bx; $bx = $cx; $cx = $u;
	$fa = $fb; $fb = $fc; $fc = $fu;
    }

    return ($ax, $bx, $cx, $fa, $fb, $fc);
}

my $CGOLD= 0.3819660;
my $ZEPS= 1e-10;

sub Brent
{
    my ($ax, $bx, $cx, $func, $tol, $ITMAX) = @_;
    my ($d, $u, $x, $w, $v); # ordinates
    my ($fu, $fx, $fw, $fv); # function evaluations

    $ITMAX = 100 unless (defined $ITMAX);
    $tol = 1e-8 unless (defined $tol);

    my $a = min($ax, $cx);
    my $b = max($ax, $cx);

    $x = $w = $v = $bx;
    $fx = $fw = $fv = &$func($x);
    my $e = 0.0; # will be distance moved on the step before last
    my $iter = 0;

    while ($iter < $ITMAX) 
    {
	my $xm = 0.5 * ($a + $b);
	my $tol1 = $tol * abs($x) + $ZEPS;
	my $tol2 = 2.0 * $tol1;

	last if (abs($x - $xm) <= ($tol2 - 0.5 * ($b - $a)));

	if (abs($e) > $tol1)
	{
	    my $r = ($x-$w)*($fx-$fv);
	    my $q = ($x-$v)*($fx-$fw);
	    my $p = ($x-$v)*$q-($x-$w)*$r;
	    if (($q = 2*($q-$r)) > 0.0) { $p = -$p; }
	    $q = abs($q);
	    my $etemp = $e;
	    $e = $d;

	    if ( (abs($p) >= abs(0.5 * $q * $etemp)) ||
		($p <= $q*($a - $x)) || ($p >= $q * ($b - $x)) )
	    {
		goto gsec;
	    }

	    # parabolic step OK here - take it
	    $d = $p/$q;
	    $u = $x + $d;

	    if ( (($u - $a) < $tol2) || (($b - $u) < $tol2) )
	    {
		$d = sign($tol1, $xm - $x);
	    }
	    goto dcomp;
	}

      gsec: # We arrive here for a Golden section step
	$e = (($x >= $xm) ? $a : $b) - $x;
	$d = $CGOLD * $e; # Golden section step

      dcomp: # We arrive here with d from Golden section or parabolic step
	$u = $x + ((abs($d) >= $tol1) ? $d : sign($tol1, $d));
	$fu = &$func($u); # 1 &$function evaluation per iteration

	#
	# Decide what to do with &$function evaluation
	#
	if ($fu <= $fx)
	{
	    if ($u >= $x)
	    {
		    $a = $x;
	    }
	    else
	    {
		    $b = $x;
	    }
	    $v = $w; $fv = $fw;
	    $w = $x; $fw = $fx;
	    $x = $u; $fx = $fu;
	}			
	else
	{
	    if ($u < $x)
	    {
		    $a = $u;
	    }
	    else
	    {
		    $b = $u;
	    }

	    if ($fu <= $fw || $w == $x)
	    {
		$v = $w; $fv = $fw;
		$w = $u; $fw = $fu;
	    }
	    elsif ( $fu <= $fv || $v == $x || $v == $w )
	    {
		    $v = $u; $fv = $fu;
	    }
	}

	$iter++;
    }

    carp "Brent Exceed Maximum Iterations.\n" if ($iter >= $ITMAX);
    return ($x, $fx);
}


1;
