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
  {
    printf("Debug information:\n");
    char* oToStr = OrePolyToString(P);
    printf("The calculated P was: %s\n", oToStr);
    free(oToStr);
    oToStr = OrePolyToString(f);
    printf("The claimed P^2 + P + 1 was: %s\n",oToStr);
    free(oToStr);
    oToStr = OrePolyToString(f_P);
    printf("The claimed product f*P was: %s\n",oToStr);
    free(oToStr);
    oToStr = OrePolyToString(P_f);
    printf("The claimed product P*f was: %s\n",oToStr);
    free(oToStr);
    testSuccess = 0;
  }
  free(f);
  free(P);
  free(f_P);
  free(P_f);
  return testSuccess;
}//test_1

/**
The following tests check the basic functionality of the Galois-field operations.
This includes:
1. Adding two elements in the GF
   1.a Adding two zeros
   1.b Adding a zero and a non-zero.
   1.c Adding two non-zero-elements.
 */
int test_GF_functionality()
{//test_GF_functionality
  int testSuccess = 1;
  char *temp;
  struct GFModulus z = getZeroElem();
  struct GFModulus k = addGF(z,z);
  //1.a
  if (isZero_GF(k)==0)
  {
    temp = GFModulusToString(k);
    printf("1.a failed. Result of adding two zeros should be zero. Result was: %s\n",temp);
    free(temp);
    testSuccess = 0;
  }
  //1.b
  struct GFModulus a;
  int i;
  for (i = 0; i<DEGREEEXTENSION; ++i)
  {
    a.coeffs[i] = 1;
  }
  k = addGF(z,a);
  if (isEqual_GF(a,k) ==0)
  {
    temp = GFModulusToString(k);
    printf("1.b failed. Result of adding a zero to a non-zero should not change anything. Result was: %s\n",temp);
    free(temp);
    testSuccess = 0;
  }
  struct GFModulus b;
  for (i = 0; i<DEGREEEXTENSION; ++i)
  {
    b.coeffs[i] = 2;
  }
  k = addGF(a,a);
  if (isEqual_GF(k,b)==0)
  {
    temp = GFModulusToString(k);
    printf("1.c failed. Result of adding two non-zero polynomials was wrong. Result was: %s\n", temp);
    free(temp);
    testSuccess = 0;
  }
  return testSuccess;
}//test_GF_functionality

int main()
{//main
//  int i;
  srand(time(NULL));
  if (test_1() == 0)
    printf("test_1 failed.\n");
  else
    printf("test_1 succeeded.\n");
  if (test_GF_functionality() == 0)
    printf("test_GF_functionality failed.\n");
  else
    printf("test_GF_functionality succeeded.\n");
  return 0;
}//main
