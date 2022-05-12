###############################################################################
# Latex Thesis Builder                                                        #
# Einar Rasmussen
#
# Get started: `make help`
#
# Assumes you add files you care about to a git repo.
# Need pacakages:
# ....
# pacaur -S opendetex
# Developed against GNU Make
# Last updated: 2021-06-02
###############################################################################


#########
# FLAGS #
#########

CC=lualatex
BB=biber
DETEX=detex #opendetex #detex
OUTPUT=master.pdf
FINAL=final.pdf
DRAFT=draft.pdf
VOCAB=vocabulary.dict
ASVOCAB=aspell_vocabulary.dict
CFLAGS=--output-directory=build -file-line-error -halt-on-error -interaction=nonstopmode -shell-escape
BFLAGS=--output-directory=build #--output_safechars #--output-encoding=UTF-8
NOSPELLCHECK:=report.tex preamble.tex minted.tex collaboration.tex files.tex draft.tex final.tex preamble.tex master.tex
USERTEX:=$(sort $(filter-out $(NOSPELLCHECK), $(wildcard sections/*.tex)))
DATE:=$(shell date -u +"%Y-%m-%dT%H:%M:%SZ")
# See write-good options here # https://github.com/btford/write-good#checks
# The format for the --add-tex-command is documented here:
# http://aspell.net/man-html/Notes-on-Various-Filters-and-Filter-Modes.html#TeX_002fLaTeX-Filter
ASPELLFLAGS:=--mode=tex --lang=en_US --tex-check-comments --home-dir=. --personal=$(ASVOCAB) --add-tex-command="eqref p:texttt p:subsubsection p"
ISPELLFLAGS:=-t -p $(VOCAB)
HUNSPELLFLAGS:=-t -p $(VOCAB) #-r
WRITEGOODFLAGS:=--no-passive #--yes-eprime
LANGUAGETOOLFLAGS:=--disable WHITESPACE_RULE

###########
# Aliases #
###########

# first rule in file means it is executed when you just write "make"
default: final

# makefile-rules that don't produce files directly
.PHONY: pdf zip xz usage spell clean texcount wc draft final

help: usage

d: draft

#########
# Rules #
#########

# Prints usage message
usage:
	@echo 'The following targets are available:'
	@echo 'make            compile the final document with bibliography.'
	@echo 'make draft      compile a draft document as fast as possible.'
	@echo 'make clean      clean out the build/ directory.  Uses git -X.'
	@echo 'make texcount   count the lines, words and characters in the document.'
	@echo 'make proofing   runs all proofing tools below.'
	@echo 'make spell      runs spellchecker.'
	@echo 'make grammar    runs write-good grammar checker.'
	@echo 'make chktex     runs chktex tex checker.'
	@echo 'make lix        calculate lix.'
	@echo 'make zip        collect all files in a zip archive.  Uses git archive to decide inclusion.'

# This recipe is documented here
# https://en.wikibooks.org/wiki/LaTeX/Bibliography_Management#Why_won't_LaTeX_generate_any_output?
final: *.tex *.bib Makefile #clean
	-rm build/$(OUTPUT:.pdf=.tex)
	ln -sf final.tex master.tex
	$(CC) $(CFLAGS) $(OUTPUT:.pdf=.tex)
	$(BB) $(BFLAGS) $(basename $(OUTPUT))
	$(CC) $(CFLAGS) $(OUTPUT:.pdf=.tex)
	$(CC) $(CFLAGS) $(OUTPUT:.pdf=.tex)
	#@clear
	texcount $(USERTEX)
	@echo "Compilation succeeded."

# Build faster

draft: *.tex
	ln -sf draft.tex $(OUTPUT:.pdf=.tex)
	$(CC) $(CFLAGS)  $(OUTPUT:.pdf=.tex)
	#clear
	texcount $(USERTEX)
	@echo "Compilation succeeded."

# Clean based on git
clean:
	-mv missfont.log build/
	git clean -X --interactive build/ #  --dry-run
	@#rm --verbose !(.gitignore) build/*  Does not work.
	@#latexrun --clean-all

############
# Proofing #
############

proofing: chktex spell grammar

# Select best available spell- and grammarchecker.
spell: hunspell
grammar: languagetool

chktex: $(USERTEX)
	chktex $^

# interactive
ispell: $(USERTEX)
	@$(foreach FILE, $^, \
	    echo "-----" $(FILE) "-----"; \
	    ispell -p $(ISPELLFLAGS) $(FILE);)
	@echo "Please correct or append these words to your user dictionary."
	@echo "Spellcheck with ispell completed."

ispelllist: $(USERTEX)
	@$(foreach FILE, $^, \
	    echo "-----" $(FILE) "-----"; \
	    $(DETEX) $(FILE) | ispell $(ISPELLFLAGS) -l;)
	@echo "Please correct or append these words to your user dictionary $(VOCAB)."
	@echo "Spellcheck with ispell completed."

aspell: $(USERTEX)
	@$(foreach FILE, $^, \
	    echo "-----" $(FILE) "-----"; \
	    aspell $(ASPELLFLAGS) check $(FILE);)
	@echo "Spellcheck with aspell completed."

aspelllist: $(USERTEX)
	@$(foreach FILE, $^, \
	    echo "-----" $(FILE) "-----"; \
	    $(DETEX) $(FILE) | aspell $(ASPELLFLAGS) list;)
	@echo "Please correct or append these words to your user $(VOCAB)."
	@echo "Spellcheck with aspell completed."

hunspell: $(USERTEX)
	@$(foreach FILE, $^, \
	    echo "-----" $(FILE) "-----"; \
	    hunspell $(HUNSPELLFLAGS) $(FILE);)
	@echo "Spellcheck with hunspell completed."

hunspelllist: $(USERTEX)
	@$(foreach FILE, $^, \
	    echo "-----" $(FILE) "-----"; \
	    hunspell $(HUNSPELLFLAGS) -l $(FILE);)
	@echo "Please correct or append these words to your user dictionary $(VOCAB)."
	@echo "Spellcheck with hunspell completed."

write-good: $(USERTEX)
	@$(foreach FILE, $^, \
	    echo "-----" $(FILE) "-----"; \
	    $(DETEX) $(FILE) > build/write-good.temp; \
	    npx write-good $(WRITEGOODFLAGS) build/write-good.temp; \
	    rm build/write-good.temp;)

pylanguagetool: $(USERTEX)
	@echo "PRIVACY WARNING: This sends your data to 3rd party.  To continue, run 'make languagetool consent=yes'."
	test $(consent)
	$(foreach FILE, $^, \
	    echo "\n-----" $(FILE) "-----\n"; \
	    $(DETEX) $(FILE) | pylanguagetool ;)

languagetool: $(USERTEX)
	@$(foreach FILE, $^, \
	    echo "\n-----" $(FILE) "-----\n"; \
	    $(DETEX) $(FILE) > build/languagetool.temp; \
	    languagetool $(LANGUAGETOOLFLAGS) build/languagetool.temp; \
	    rm build/languagetool.temp;)

##############
# Statistics #
##############

texcount: $(USERTEX)
	texcount $^

wc: $(USERTEX)
	cat $^ | $(DETEX) | wc

lix:
	echo "Not implemented"

# for i in {1..3}; do cat $i.tex | $(DETEX) | tr -d '[:space:]' |  wc --chars; done
# cat ?.tex | $(DETEX) | tr -d '[:space:]' |  wc --chars | tr -d '\n' | cat <(echo -n "scale=2;") - <(echo -n "/(2400)\n")  | tee /dev/tty | bc -ql
#| cat `echo -n "scale=2;"` - `echo -n "/2400\n"`
#| tee /dev/tty
normalsider: ?.tex
	echo "Fordeling af normalsider à 2400 anslag uden mellemrum:"
	@$(foreach FILE, $(sort $^), \
	    echo -n  $(FILE) ": "; \
	    $(DETEX) $(FILE) \
	    | tr -d '[:space:]' \
	    | wc --chars \
	    | tr -d '\n' \
	    | awk '{print "scale=2;" $$1 "/2400"}' \
	    | bc -ql;)
	$(MAKE) totalnormalsider

totalnormalsider: ?.tex
	echo -n "Total : "
	$(DETEX) $^ \
	| tr -d '[:space:]' \
	| wc --chars \
	| tr -d '\n' \
	| awk '{print "scale=2;" $$1 "/2400"}' \
	| bc -ql

#Tæl normalside: 2400 ord/s ifølge Pernille Bjørns UIS2018 slides til første forelæsnig.
uispages: ?.tex
	echo "Fordeling af normalsider à 400 ord:"
	@$(foreach FILE, $(sort $^), \
	    echo -n  $(FILE) ": "; \
	    $(DETEX) $(FILE) \
	    | wc --words \
	    | tr -d '\n' \
	    | awk '{print "scale=2;" $$1 "/400"}' \
	    | bc -ql;)


#############
# Packaging #
#############

zip: *.tex *.bib Makefile $(OUTPUT)
	git archive --format=zip -9 --verbose --output=exports/$(DATE).thesis.zip HEAD
	# Add print
	zip -j exports/$(DATE).thesis.zip build/master.pdf
	# Test integrity
	zip --test exports/$(DATE).thesis.zip

xz: zip
	# extract zip and pipe into xz
	unzip -p exports/$(DATE).thesis.zip \
	| xz --compress -9 --extreme --verbose \
	> exports/$(DATE).thesis.xz
	# Test integrity
	xz --test exports/$(DATE).thesis.xz


################
# Experimental #
################

# Validate the bibtex bibliography
validatebib: build/master.bcf
	biber --validate_datamodel $?

# Remove metadata from PDF
anonymize:
	#https://gist.github.com/hubgit/6078384

latexmk:
	echo "Not implemented"

latexrun:
	echo "Not implemented"

progress: $(USERTEX)
	cat $^ | $(DETEX) | tr -d | wc --chars | cat <(echo -n " 100/( 4 * 2400 ) * ") - | tee /dev/tty | bc -ql
