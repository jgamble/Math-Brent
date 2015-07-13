use 5.008003;
use strict;
use warnings;
use Module::Build;

my $builder = Module::Build->new(
	module_name   => 'Math::Brent',

	dist_abstract => 'Single Dimensional Function Minimisation',
	dist_author   => [q{J. A. R. Williams},
			q{John M. Gamble <jgamble@cpan.org>}],
	dist_version => '0.04',
	dist_name => 'Math-Brent',

	requires => {
		perl => '5.8.3',
		'Math::Fortran' => 0,
		'Math::VecStat' => 0
	},

	configure_requires => {
		'Module::Build' => 0,
	},

	build_requires => {
		'Test::More' => 0,
	},

	license => 'perl',
	create_license => 1,
	create_readme => 0,
	create_makefile_pl => 'traditional',
	dynamic_config => 0,
	release_status => 'stable',
	add_to_cleanup => [ 'Math-Brent-*' ],
);

$builder->create_build_script();
