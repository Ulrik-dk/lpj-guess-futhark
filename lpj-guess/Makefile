all: 0 1 2 3

.PHONY: 3

r: #results
	gprof example/guess example/gmon.out > example/analysis.txt
	less example/analysis.txt

0:
	futhark check futsource/framework/guess_datatypes.fut
	futhark check futsource/framework/guess.fut
	futhark check futsource/framework/guessmath.fut
	futhark check futsource/framework/parameters.fut

	#futhark check futsource/modules/bvoc.fut
	futhark check futsource/modules/ncompete.fut
	futhark check futsource/modules/q10.fut
	futhark check futsource/modules/soil.fut
	futhark check futsource/modules/canexch.fut

	futhark check futsource/everything.fut
	futhark check futsource/futhark-extras.fut

1:
	sh compile.sh

2:
	rm -rf $(FUTHARK_TEST_DIR)
	mkdir $(FUTHARK_TEST_DIR)
	cd example;	sh lpjguess_run_cf.sh

FUTHARK_TEST_DIR = ./example/futhark_tests
3:
	$(foreach file,	$(wildcard	$(FUTHARK_TEST_DIR)/*.fut),	futhark test $(file);)
