#This is the makefile of the ore_algebra package
#(c) by Reinhold Burger and Albert Heinle


#Compiler
CC=gcc
CFLAGS=-Wall
CPFLAGS=-Wall -fopenmp
COPFLAGS=-Wall -fopenmp -c
COFLAGS=-Wall -c

all: impl tests

tests: ore_algebra.o gf_coefficients.o
	$(CC) $(CFLAGS) ore_algebra_test.c -o tests ore_algebra.o gf_coefficients.o

impl: ore_algebra.o gf_coefficients.o
	$(CC) $(CFLAGS) impl.c -o impl ore_algebra.o gf_coefficients.o

ore_algebra.o:
	$(CC) $(COFLAGS) ore_algebra.c -o ore_algebra.o

gf_coefficients.o:
	$(CC) $(COFLAGS) gf_coefficients.c -o gf_coefficients.o

clean:
	rm -f gf_coefficients.o
	rm -f ore_algebra.o
	rm -f impl
	rm -f tests
