=head1 NAME

Math::Brent - Single Dimensional Function Minimisation

=head1 SYNOPSIS

    use Math::Brent qw(Minimise1D);

    my ($x, $y) = Minimise1D($guess, $scale, \&func, $tol, $itmax);

or

    use Math::Brent qw(BracketMinimum Brent);

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

=head3 Minimise1D()

Provides a simple interface to the L</BracketMinimum()> and L</Brent()>
routines.

Given a function, an initial guess for the function's
minimum, and its scaling, this routine converges
to the function's minimum using Brent's method.

    ($x, $y) = Minimise1D($guess, $scale, \&func);

The minimum is reached within a certain tolerance (defaulting 1e-7), and
attempts to do so within a maximum number of iterations (defaulting to 100).
You may override them by providing alternate values:

    ($x, $y) = Minimise1D($guess, $scale, \&func, 1.5e-8, 120);

=head3 Brent()

Given a function and a triplet of abcissas B<$ax>, B<$bx>, B<$cx>, such that

=over 4

=item 1. B<$bx> is between B<$ax> and B<$cx>, and

=item 2. B<func($bx)> is less than both B<func($ax)> and B<func($cx)>),

=back

Brent() isolates the minimum to a fractional precision of about B<$tol>
using Brent's method.

A maximum number of iterations B<$itmax> may be specified for this search - it
defaults to 100. Returned is a list consisting of the abcissa of the minum
and the function value there.

=head3 BracketMinimum()

Given a function reference B<\&func> and
distinct initial points B<$ax> and B<$bx> searches in the downhill
direction (defined by the function as evaluated at the initial points)
and returns a list of the three points B<$ax>, B<$bx>, B<$cx> which
bracket the minimum of the function and the function values at those
points.

=head1 EXAMPLE

    use Math::Brent qw(Minimise1D);

    sub sinc {
      my $x = shift ;
      return $x ? sin($x)/$x: 1;
    }

    my($x, $y) = Minimise1D(1, 1, \&sinc, 1e-7);
    print "Minimum is at sinc($x) = $y\n";

produces the output

    Minimum is at sinc(4.4934094397196) = -.217233628211222

Anonymous subroutines may also be used as the function reference:

    my $cubic_ref = sub {my($x) = @_; return 6.25 + $x*$x*(-24 + $x*8));};

    my($x, $y) = Minimise1D(3, 1, $cubic_ref);
    print "Minimum of the cubic at $x = $y\n";


=head1 BUGS

Please report any bugs or feature requests via Github's L<issues link|:q>

=head1 AUTHOR

John A.R. Williams B<J.A.R.Williams@aston.ac.uk>

John M. Gamble B<jgamble@cpan.org> (current maintainer)

=head1 SEE ALSO

"Numerical Recipies: The Art of Scientific Computing"
W.H. Press, B.P. Flannery, S.A. Teukolsky, W.T. Vetterling.
Cambridge University Press. ISBN 0 521 30811 9.

Richard P. Brent, L<Algorithms for Minimization Without Derivatives|http://www.worldcat.org/title/algorithms-for-minimization-without-derivatives/oclc/515987&referer=brief_results>

Professor (Emeritus) Richard Brent has a web page at
L<http://maths-people.anu.edu.au/~brent/>

=cut

package Math::Brent;

use strict;
use warnings;
use 5.10.1;

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

our $VERSION = 0.05;

use Math::VecStat qw(max min);
use Math::Utils qw(:fortran);
use Carp;

sub Minimise1D
{
    my ($guess, $scale, $func, $tol, $itmax) = @_;
    my ($a, $b, $c) = BracketMinimum($guess - $scale, $guess + $scale, $func);

    return Brent($a, $b, $c, $func, $tol, $itmax);
}

#
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

my $GOLD = 0.5 + sqrt(1.25); # Default magnification ratio for intervals is phi.
my $GLIMIT = 100.0; # Max magnification for parabolic fit step
my $TINY = 1E-20;

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
	my $t = $ax; $ax = $bx; $bx = $t;
	$t = $fa; $fa = $fb; $fb = $t;
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
	    (2.0 * copysign(max(abs($q - $r), $TINY), $q - $r));

	my $ulim = $bx + $GLIMIT * ($cx - $bx); # We won't go further than this
	my $fu;

	#
	# Parabolic U between B & C - try it
	#
	if (($bx - $u) * ($u - $cx) > 0.0)
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
		$u = $cx + $GOLD * ($cx - $bx);
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
	    # eject parabolic U, use default magnification
	    $u = $cx + $GOLD * ($cx - $bx);
	    $fu = &$func($u);
	}

	# Eliminate oldest point & continue
	$ax = $bx; $bx = $cx; $cx = $u;
	$fa = $fb; $fb = $fc; $fc = $fu;
    }

    return ($ax, $bx, $cx, $fa, $fb, $fc);
}

#
# The complementary step is (3 - sqrt(5))/2, which resolves to 2 - phi.
#
my $CGOLD = 2 - $GOLD;
my $ZEPS = 1e-10;

sub Brent
{
    my ($ax, $bx, $cx, $func, $tol, $ITMAX) = @_;
    my ($d, $u, $x, $w, $v); # ordinates
    my ($fu, $fx, $fw, $fv); # function evaluations

    $ITMAX //= 100;
    $tol //= 1e-8;

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
	    my $r = ($x-$w) * ($fx-$fv);
	    my $q = ($x-$v) * ($fx-$fw);
	    my $p = ($x-$v) * $q-($x-$w)*$r;

	    $p = -$p if (($q = 2 * ($q - $r)) > 0.0);

	    $q = abs($q);
	    my $etemp = $e;
	    $e = $d;

	    unless ( (abs($p) >= abs(0.5 * $q * $etemp)) ||
		($p <= $q * ($a - $x)) || ($p >= $q * ($b - $x)) )
	    {
                #
	        # Parabolic step OK here - take it.
                #
	        $d = $p/$q;
	        $u = $x + $d;

	        if ( (($u - $a) < $tol2) || (($b - $u) < $tol2) )
	        {
		    $d = copysign($tol1, $xm - $x);
	        }
	        goto dcomp; # Skip the golden section step.
	    }
	}

        #
        # Golden section step.
        #
	$e = (($x >= $xm) ? $a : $b) - $x;
	$d = $CGOLD * $e;

        #
        # We arrive here with d from Golden section or parabolic step.
        #
        dcomp:
	$u = $x + ((abs($d) >= $tol1) ? $d : copysign($tol1, $d));
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
