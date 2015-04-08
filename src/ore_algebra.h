/*
  Header file for ore_algebra.c

  (c) Reinhold Burger and Albert Heinle
 */

#ifndef ORE_ALGEBRA_H
#define ORE_ALGEBRA_H

#include "gf_coefficients.h"

#define MIN_TOTAL_DEG_FOR_PARALLEL (200)

struct OrePoly
{//OrePoly
  int degD1;//representing the maximal degree in D1 appearing in the
	    //given Ore Polynomial
  int degD2;//representing the maximal degree in D2 appearing in the
	    //given Ore Polynomial
  struct GFModulus* coeffs;//array of size (degD1+1)*(degD2+1)
  struct GFModulus (*ptrD1manip)(struct GFModulus); //Sets the commutation rules for the ground
			  //field with d1
  struct GFModulus (*ptrD2manip)(struct GFModulus); //Sets the commutation rules for the ground
			  //field with d2
};//OrePoly

char* OrePolyToString(struct OrePoly*);

void OrePolyToStdOut(struct OrePoly*);

struct OrePoly * getOrePolyViaIntegerCoefficients(int, int,struct GFModulus (*ptrD1manip)(struct GFModulus),
                                                  struct GFModulus(*ptrD2manip)(struct GFModulus),int*);


int isZero_OrePoly(struct OrePoly*);

int isEqual_OrePoly(struct OrePoly*, struct OrePoly*);

struct OrePoly *getIdentityElemOrePoly(void);

struct OrePoly *getZeroElemOrePoly(void);

struct OrePoly* scalarMult(int, struct OrePoly*);

struct OrePoly* add(struct OrePoly*, struct OrePoly*);

struct OrePoly* minus(struct OrePoly*, struct OrePoly*);

struct OrePoly* mult(struct OrePoly*, struct OrePoly*);

struct OrePoly * getRandomOrePoly(int, int, struct GFModulus(*ptrD1manip)(struct GFModulus),
				  struct GFModulus (*ptrD2manip)(struct GFModulus));

struct OrePoly * generateRandomSecretKey(int, struct OrePoly*);


#endif
