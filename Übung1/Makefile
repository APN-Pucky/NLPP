# program name
PROGRAM = main.out
PWD = $(shell pwd)
ZIPFILE = $(shell basename $(PWD))
DATE = $(shell date)

# compiler
FC := gfortran

# compile
FCFLAGS = -Wall -Wextra -g -fbacktrace -pedantic -Wconversion #-x f95 -std=f2008
# link
FLFLAGS =

# source files and objects
#reverse = $(if $(wordlist 2,2,$(1)),$(call reverse,$(wordlist 2,$(words $(1)),$(1))) $(firstword $(1)),$(1))
SRC_O = $(patsubst %.f90, %.o, $(wildcard *.f90))
SRC_MOD = $(patsubst %.f90, %.o, $(wildcard module/*.f90))
SRC_OUT = $(patsubst %.f90, %.out, $(wildcard *.f90))
TXTS = $(patsubst %.txt, %.pdf, $(wildcard *.txt))

default: one

one: $(PROGRAM)

$(PROGRAM): $(SRC_MOD) $(SRC_O) # make single.out from all .f90
	$(FC) $(FLFLAGS) -o $@ $^

all: $(SRC_OUT) #make single .f90 to single .out

test: one
	./$(PROGRAM)


zip: clean pdfs
	cd ..; zip $(PWD)/$(ZIPFILE) -r $(ZIPFILE)
	


pdfs: $(TXTS)


%.pdf: %.txt
	sed -i 's/{dir}/$(ZIPFILE)/g' *.txt
	sed -i 's/{date}/$(DATE)/g' *.txt
	libreoffice --convert-to "pdf" $<
	#text2pdf $< #name .pdf => .txt.pdf

%.o: %.f90
	$(FC) $(FCFLAGS) -o $@ -c $<

%.out: %.f90
	$(FC) $(FCFLAGS) -o $@ $<


#%.mo: %.h
#	$(FC) $(FCFLAGS) -o $@ -c $<

clean:
	rm -f **/*.o *.o **/*.mod *.mod $(PROGRAM) $(ZIPFILE).zip
