all: build run-tests r

.PHONY: 3

build:
	sh compile.sh

run-tests:
	cd example;	sh lpjguess_run_cf.sh

r: #results
	gprof example/guess example/gmon.out > example/analysis.txt
	less example/analysis.txt

archdependencies:
	yay -S cmake netcdf

0:
	futhark check futsource/futhark-extras.fut

	futhark check futsource/framework/guessmath.fut
	futhark check futsource/framework/parameters.fut
	futhark check futsource/framework/guess_datatypes.fut
	futhark check futsource/framework/guess.fut


	futhark check futsource/modules/q10.fut
	futhark check futsource/modules/soil.fut
	futhark check futsource/modules/canexch.fut
	
	futhark check futsource/everything.fut

1: build

2:
	rm -rf $(FUTHARK_TEST_DIR)
	mkdir $(FUTHARK_TEST_DIR)
	make run-tests

FUTHARK_TEST_DIR = ./example/futhark_tests
3:
	$(foreach file,	$(wildcard	$(FUTHARK_TEST_DIR)/*.fut),	futhark test $(file);)
