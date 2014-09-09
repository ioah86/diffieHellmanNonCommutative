/*
  Test of the implementation of the Diffie-Hellman protocol -- Non-commutative version

  (c) Reinhold Burger and Albert Heinle
 */

#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include "ore_algebra.c"

/**
   This test covers the following scenario:
   Generate a random Ore-Polynomial P, then calculate
   f(P) for f = x^2 + x +1. After that, check if f(P)*P == P*f(P)
 */
int test_1()
{//test_1
  struct OrePoly* P =
    getRandomOrePoly(10,10,&Hom1,&Hom2);
  struct OrePoly *f = malloc(sizeof(struct OrePoly));
  f->degD1 = 0;
  f->degD2 = 0;
  f->ptrD1manip = &Hom1;
  f->ptrD2manip = &Hom2;
  struct GFModulus coeffsForF[1];
  coeffsForF[0].coeffs[0] = 1;
  coeffsForF[0].coeffs[1] = 0;
  coeffsForF[0].coeffs[2] = 0;
  f->coeffs = coeffsForF;
  struct OrePoly * tempf;
  struct OrePoly *tempP;
  tempf = add(f,P);
  free(f);
  f = tempf;
  tempP = mult(P,P);
  tempf = add (f, tempP);
  free(f);
  f = tempf;
  free(tempP);

  //Now, f(P) and P are generated.
  struct OrePoly *f_P = mult(f,P);
  struct OrePoly *P_f = mult(P,f);
  
  int testSuccess = 1;
  if (isEqual_OrePoly(*f_P, *P_f)==0)
    testSuccess = 0;
  free(f);
  free(P);
  free(f_P);
  free(P_f);
  return testSuccess;
}//text_1

int main()
{//main
//  int i;
  srand(time(NULL));
  if (test_1() == 0)
    printf("test_1 failed.\n");
  else
    printf("test_1 succeeded.\n");
  
  return 0;
}//main
