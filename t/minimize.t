use Test::More tests => 1;

use Math::Brent qw(:all);
use strict;
use warnings;

my $cubic_ref = sub {my($x) = @_; return 6.25 + $x*$x*(-24 + $x*8);};

my ($x, $y) = Minimise1D(3, 0.5, $cubic_ref);

ok((($y + 0.005 < 2.907589) && ($y + 0.005 < 2.907589)), "Cubic test");

1;
