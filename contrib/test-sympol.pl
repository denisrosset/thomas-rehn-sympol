#!/usr/bin/perl

use strict;
use warnings;
use File::Find;

my $sympol_root = '../build/debug';

if ($#ARGV >= 0) {
    $sympol_root = '../build/' . $ARGV[0];
}
print "using root $sympol_root\n";

sub test_polyhedron {
	my ($parameters, $filename, $expected_rays) = @_;
	print "$filename \@ $parameters";
	my $cmd = "${sympol_root}/sympol/sympol -t $parameters -i ${sympol_root}/data/$filename";
    my $sympol = `$cmd`;
	my $time = '';
	if ($sympol =~ /elapsed time: (\d+\.?\d*) seconds/) {
		$time = $1;
	}
	if ($sympol =~ /(\d+)\s+(\d+)\s+rational/) {
		my $number_of_rays = $1;
		my $dimension = $2;
		if ($number_of_rays != $expected_rays) {
			warn "  failed: $number_of_rays != $expected_rays";
		} else {
			print "  PASS  $time";
		}
	} else {
		warn "  output format mismatch";
		print $sympol;
	}
	print "\n";
}

test_polyhedron('-d', 'misc/santos_prismatoid.ext', 6);
test_polyhedron('-d --cdd', 'misc/santos_prismatoid.ext', 6);
test_polyhedron('-d', 'misc/santos_prismatoid-reduced_symmetry.ext', 12);
test_polyhedron('-d --cdd', 'misc/santos_prismatoid-reduced_symmetry.ext', 12);

test_polyhedron('-d --cdd', 'metric/metric_5.ine', 2);
test_polyhedron('-d', 'metric/metric_5.ine', 2);

test_polyhedron('-d', 'cyclic/cyclic4-5.ext', 1);
test_polyhedron('-d --cdd', 'cyclic/cyclic4-5.ext', 1);
# homogenized polar .ine contains cone apex
test_polyhedron('-d', 'cyclic/cyclic4-5.ine', 2);

test_polyhedron('-d --cdd', 'metric/metric_6.ine', 3);
test_polyhedron('-a 10 --cdd', 'metric/metric_6.ine', 3);
test_polyhedron('-a 10', 'metric/metric_6.ine', 3);
test_polyhedron('--idm-adm 5 10', 'metric/metric_6.ine', 3);
test_polyhedron('--idm-adm-level 1 2', 'metric/metric_6.ine', 3);
test_polyhedron('--cdd --idm-adm-level 1 1', 'metric/metric_6.ine', 3);

test_polyhedron('-d --cdd', 'voronoi_cones/d4.ine', 2);
test_polyhedron('-d', 'voronoi_cones/d4.ine', 2);

test_polyhedron('-d --cdd', 'voronoi_cones/d5.ine', 3);
test_polyhedron('-d', 'voronoi_cones/d5.ine', 3);

test_polyhedron('-d', 'voronoi_cones/e6.ine', 12);
test_polyhedron('-d --cdd', 'voronoi_cones/e6.ine', 12);
test_polyhedron('-a 20', 'voronoi_cones/e6.ine', 12);
test_polyhedron('--idm-adm 5 10', 'voronoi_cones/e6.ine', 12);
test_polyhedron('--idm-adm-level 0 2', 'voronoi_cones/e6.ine', 12);
test_polyhedron('--idm-adm-level 0 1', 'voronoi_cones/e6.ine', 12);
test_polyhedron('--cdd --idm-adm-level 0 1', 'voronoi_cones/e6.ine', 12);
test_polyhedron('--idm-adm-level 1 2', 'voronoi_cones/e6.ine', 12);

