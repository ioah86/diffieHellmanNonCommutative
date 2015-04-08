#This is the makefile of the ore_algebra package
#(c) by Reinhold Burger and Albert Heinle


#Compiler
CC=gcc-mp-4.9
CFLAGS=-Wall
CPFLAGS=-Wall -fopenmp
COPFLAGS=-Wall -fopenmp -c
COFLAGS=-Wall -c
SRCFLDR=src
BLDFLDR=build
LIBFLDR=lib
TSTSFLDR=tests

all: impl tests

tests: ore_algebra gf_coefficients $(TSTSFLDR)/ore_algebra_test.c
	$(CC) $(CPFLAGS) $(TSTSFLDR)/ore_algebra_test.c -o $(BLDFLDR)/tests $(LIBFLDR)/ore_algebra.o $(LIBFLDR)/gf_coefficients.o

impl: ore_algebra gf_coefficients $(SRCFLDR)/impl.c
	$(CC) $(CPFLAGS) $(SRCFLDR)/impl.c -o $(BLDFLDR)/impl $(LIBFLDR)/ore_algebra.o $(LIBFLDR)/gf_coefficients.o

ore_algebra: $(SRCFLDR)/ore_algebra.c
	$(CC) $(COPFLAGS) $(SRCFLDR)/ore_algebra.c -o $(LIBFLDR)/ore_algebra.o

gf_coefficients: $(SRCFLDR)/gf_coefficients.c
	$(CC) $(COFLAGS) $(SRCFLDR)/gf_coefficients.c -o $(LIBFLDR)/gf_coefficients.o

clean:
	rm -f $(LIBFLDR)/ore_algebra.o
	rm -f $(LIBFLDR)/gf_coefficients.o
	rm -f $(BLDFLDR)/tests
	rm -f $(BLDFLDR)/impl
