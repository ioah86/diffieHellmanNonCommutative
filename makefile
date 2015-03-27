#This is the makefile of the ore_algebra package
#(c) by Reinhold Burger and Albert Heinle


#Compiler
CC=gcc
CFLAGS=-Wall
CPFLAGS=-Wall -fopenmp
COPFLAGS=-Wall -fopenmp -c
COFLAGS=-Wall -c
SRCFLDR=src
BLDFLDR=build
LIBFLDR=lib
TSTSFLDR=tests

all: buildfolder libfolder impl tests

buildfolder:
	mkdir -p $(BLDFLDR)

libfolder:
	mkdir -p $(LIBFLDR)

tests: ore_algebra gf_coefficients $(TSTSFLDR)/ore_algebra_test.c
	$(CC) $(CPFLAGS) $(TSTSFLDR)/ore_algebra_test.c -o $(BLDFLDR)/tests $(LIBFLDR)/ore_algebra.o $(LIBFLDR)/gf_coefficients.o

impl: ore_algebra gf_coefficients $(SRCFLDR)/impl.c
	$(CC) $(CPFLAGS) $(SRCFLDR)/impl.c -o $(BLDFLDR)/impl $(LIBFLDR)/ore_algebra.o $(LIBFLDR)/gf_coefficients.o

ore_algebra: $(SRCFLDR)/ore_algebra.c
	$(CC) $(COPFLAGS) $(SRCFLDR)/ore_algebra.c -o $(LIBFLDR)/ore_algebra.o

gf_coefficients: $(SRCFLDR)/gf_coefficients.c
	$(CC) $(COFLAGS) $(SRCFLDR)/gf_coefficients.c -o $(LIBFLDR)/gf_coefficients.o

clean:
	rm -rf $(LIBFLDR)
	rm -rf $(BLDFLDR)
