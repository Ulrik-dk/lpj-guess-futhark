default:

all: build setup-tests run-tests r

run-tests:
	cd example;	bash lpjguess_run_cf.sh


1: build

2: run-tests

FUTHARK_TEST_DIR = ./example/futhark_tests
3:
	$(foreach file,	$(wildcard	$(FUTHARK_TEST_DIR)/*.fut),	futhark test $(file);)

check:
	futhark check futsource/futhark-extras.fut

	futhark check futsource/framework/guess.fut
	futhark check futsource/framework/guessmath.fut
	futhark check futsource/framework/parameters.fut

	futhark check futsource/modules/q10.fut
	futhark check futsource/modules/soil.fut
	futhark check futsource/modules/canexch.fut

r: #results
	gprof example/guess example/gmon.out > example/analysis.txt
	atom example/analysis.txt

setup-tests:
	mv example/gridlist_cf.txt example/gridlist_cf.txt.BAK
	echo "26 32" > example/gridlist_cf.txt

build:
	bash compile.sh

b:
	bash compile.sh

archdependencies:
	yay -S cmake netcdf

permissions:
	chmod -R 777 .
