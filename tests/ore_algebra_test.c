/*
  Test of the implementation of the Diffie-Hellman protocol -- Non-commutative version

  (c) Reinhold Burger and Albert Heinle
 */

#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include "../src/ore_algebra.h"

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
  if (isEqual_OrePoly(f_P, P_f)==0)
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
2. Multiplying two elements in the GF
   2.a Multiplying 1 to an non-zero and non-one element
   2.b Multiplying 1 to 1
   2.c Multiplying zero to a non-zero element
   2.d Multiplying 0 to 0
   2.e Multiplying different non-unit non-zero elements.
3. Testing the isZero-Test
4. Testing the isEqual-Function
5. Minussing two elements
   5.a Minussing the same element with oneself
   5.b Minussing different elements
6. Scalar-multiplying in GF
   6.a Scalar multiplication with zero-element
   6.b Scalar multiplication with negative element
   6.c Scalar multiplication with positive element
 */
int test_GF_functionality()
{//test_GF_functionality
  int testSuccess = 1;
  char *temp;
  struct GFModulus z = getZeroElemGF();
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
  //1.c
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
  
  //2
  b.coeffs[0] =1;
  for (i = 1; i<DEGREEEXTENSION; ++i)
  {//making b the unit-element
    b.coeffs[i] = 0;
  }//making b the unit-element
  
  //random elements in the ground field
  struct GFModulus f1,f2,f3,f4,f5;
  f1.coeffs[0] = 1; f1.coeffs[1] = 2; f1.coeffs[2] = 1;
  f2.coeffs[0] = 1; f2.coeffs[1] = 3; f2.coeffs[2] = 1;
  f3.coeffs[0] = 1; f3.coeffs[1] = 3; f3.coeffs[2] = 3;
  f4.coeffs[0] = 1; f4.coeffs[1] = 0; f4.coeffs[2] = 1;
  f5.coeffs[0] = 0; f5.coeffs[1] = 3; f5.coeffs[2] = 4;
  //2.a
  k = multGF(b,f1);
  if (!isEqual_GF(f1,k))
  {
    testSuccess = 0;
    temp = GFModulusToString(k);
    printf("2.a failed. The unit element multiplied by a^2 + 2a + 1 resulted\
 unexpectedly in: %s\n",temp);
    free(temp);
  }
  //2.b
  k = multGF(b,b);
  if (!isEqual_GF(k,b))
  {
    testSuccess = 0;
    temp = GFModulusToString(k);
    printf("2.b failed The unit element multiplied by itself resulted\
 unexpectedly in: %s\n",temp);
    free(temp);
  }

  //2.c
  k = multGF (f2,z);
  if (!isEqual_GF(z,k))
  {
    testSuccess = 0;
    temp = GFModulusToString(k);
    printf("2.c failed. The zero element multiplied by a^2 + 3a + 1 resulted\
 unexpectedly in: %s\n",temp);
    free(temp);
  }

  //2.d
  k = multGF(z,z);
  if (!isEqual_GF(z,k))
  {
    testSuccess = 0;
    temp = GFModulusToString(k);
    printf("2.d failed. The zero element multiplied by itself resulted\
 unexpectedly in: %s\n",temp);
    free(temp);
  }
  //2.e
  //f4 * f5 = 2a^2 + 2a +1
  k = multGF(f4,f5);
  if (k.coeffs[0] != 1 || k.coeffs[1] != 2 || k.coeffs[2]!=2)
  {
    testSuccess = 0;
    temp = GFModulusToString(k);
    printf("2.e failed. The f4*f5 resulted unexpectedly in: %s\n",temp);
    free(temp);
  }
  //f3*f2 = 4*a^2 + a
  k = multGF(f3,f2);
  if (k.coeffs[0] != 0 || k.coeffs[1] != 1 || k.coeffs[2]!=4)
  {
    testSuccess = 0;
    temp = GFModulusToString(k);
    printf("2.e failed. The f3*f2 resulted unexpectedly in: %s\n",temp);
    free(temp);
  }
  //f1*f5 = 3*a^2 + 3*a + 2
  k = multGF(f1,f5);
  if (k.coeffs[0] != 2 || k.coeffs[1] != 3 || k.coeffs[2]!=3)
  {
    testSuccess = 0;
    temp = GFModulusToString(k);
    printf("2.e failed. The f1*f5 resulted unexpectedly in: %s\n",temp);
    free(temp);
  }
  //f4 * f2 = 4*a^2 + a + 2
  k = multGF(f4,f2);
  if (k.coeffs[0] != 2 || k.coeffs[1] != 1 || k.coeffs[2]!=4)
  {
    testSuccess = 0;
    temp = GFModulusToString(k);
    printf("2.e failed. The f4*f2 resulted unexpectedly in: %s\n",temp);
    free(temp);
  }
  //3.
  k = minusGF(f4,f4);
  if (!isZero_GF(k))
  {
    testSuccess = 0;
    temp = GFModulusToString(k);
    printf("3. failed. The operation f4-f4 resulted unexpectedly in: %s\n",temp);
    free(temp);
  }
  k = getIdentityElemGF();
  if (isZero_GF(k))
  {
    testSuccess = 0;
    temp = GFModulusToString(k);
    printf("3. failed. Unexpectedly, the identity-element %s was identified to \
be zero\n",temp);
    free(temp);
  }
  //4.
  if (!isEqual_GF(f2,f2))
  {
    testSuccess = 0;
    printf("4. failed. f2 is not equal to f2 according to the isEqual function.\n");
    free(temp);
  }
  if (isEqual_GF(f2,f1))
  {
    testSuccess = 0;
    printf("4. failed. According to the isEqual_GF function, f2 is equal to f1.\n");
    free(temp);
  }
  //5a
  k = minusGF (f1,f1);
  if (!isZero_GF(k))
  {
    testSuccess = 0;
    temp = GFModulusToString(k);
    printf("5.a failed. The f1-f1 resulted unexpectedly in: %s\n",temp);
    free(temp);
  }
  //5.b
  k = minusGF (f1,f3);
  if (k.coeffs[0] != 0 || k.coeffs[1] != 4 || k.coeffs[2] != 3)
  {
    testSuccess = 0;
    temp = GFModulusToString(k);
    printf("5.b failed. The f1-f3 resulted unexpectedly in: %s\n",temp);
    free(temp);
  }
  //6
  //6.a
  k = scalarMultGF(0,f1);
  if (!isZero_GF(k))
  {
    testSuccess = 0;
    temp = GFModulusToString(k);
    printf("6.a failed. The scalar-multiplication of f1 with zero unexpectedly \
resulted in: %s\n",temp);
    free(temp);
  }
  //6.b
  k = scalarMultGF(-1,f1);
  if (k.coeffs[0] != MODULUS-1 || k.coeffs[1] != MODULUS - 2 ||
      k.coeffs[2] != MODULUS-1)
  {
    testSuccess = 0;
    temp = GFModulusToString(k);
    printf("6.b failed. The scalar-multiplication of f1 with -1 unexpectedly \
resulted in: %s\n",temp);
    free(temp);
  }
  //6.c
  k = scalarMultGF(2,f2);
  if (k.coeffs[0] != 2 || k.coeffs[1] != 1 ||
      k.coeffs[2] != 2)
  {
    testSuccess = 0;
    temp = GFModulusToString(k);
    printf("6.c failed. The scalar-multiplication of f2 with 2 unexpectedly \
resulted in: %s\n",temp);
    free(temp);
  }
  return testSuccess;
}//test_GF_functionality

/**
 The tests in this function cover the basic arithmetics for the
 OrePoly structure. The tests include:
 1. test isZero
   1.a isZero to a zero Ore polynomial
   1.b isZero to a nonZero Ore polynomial
 2. test isEqual
   2.a compare two equal polynomials
   2.b compare two non-equal polynomials
 3. test scalar-multiplication with integer
   3.a multiply with zero-element
   3.b multiply with 1
   3.c multiply with non-unit non-zero element
 4. test addition
   4.a Add zero-element to itself
   4.b Add zero-element to non-zero element
   4.c Add two non-zero elements
 5. testing the multiplication
   5.a Multiply zero to a non-zero OrePoly
   5.b Multiply zero to zero
   5.c Multiply one to itself
   5.d Multiply one to a non-one non-unit element
   5.e Multiply two non-zero non-unit elements
 */
int test_OrePoly_Arith()
{//test_OrePoly_Arith
  int testSuccess = 1;
  char *tempOutp;
  struct OrePoly *z = getZeroElemOrePoly();
  struct OrePoly *one = getIdentityElemOrePoly();
  //In what follows, we will generate some random elements for
  //testing purposes
  //f1 =(a^2+a+3)*d1^2+(4*a^2+4)*d1*d2+(a^2+4*a)*d2^2+(3*a^2+4*a+4)*d1+(a^2+3)*d2
  int tempCoeffsf1[9*DEGREEEXTENSION] = {0,0,0, 4,4,3, 3,1,1, 3,0,1,
					 4,0,4, 0,0,0, 0,4,1, 0,0,0,
					 0,0,0};
  struct OrePoly *f1 =
    getOrePolyViaIntegerCoefficients(2,2,&Hom1,&Hom2,tempCoeffsf1);
  // f2=(-2*a^2-2*a-2)*d1^2+(-2*a^2+a-2)*d1*d2+(-a+1)*d2^2+(2*a^2+2*a-2)*d1+(-2*a+1)*d2
  int tempCoeffsf2[9*DEGREEEXTENSION] = {0,0,0, -2,2,2, -2,-2,-2,
					 1,-2,0, -2,1,-2, 0,0,0, 
					 1,-1,0, 0,0,0, 0,0,0};
  struct OrePoly *f2 =
    getOrePolyViaIntegerCoefficients(2,2,&Hom1,&Hom2, tempCoeffsf2);
  // f3=(2*a^2-2*a-1)*d1^2+(a^2+2)*d2^2+(a^2+2*a)*d1+(-2*a^2-2)*d2+1
  int tempCoeffsf3[9*DEGREEEXTENSION] = {1,0,0, 0,2,1, -1,-2,2,
					 -2,0,-2, 0,0,0, 0,0,0, 
					 2,0,1, 0,0,0, 0,0,0};
  struct OrePoly *f3 =
    getOrePolyViaIntegerCoefficients(2,2,&Hom1,&Hom2, tempCoeffsf3);

  //f1*f2 = (2*a^2+3*a)*d1^4+(4*a^2+4*a+3)*d1^3*d2+(3*a^2+4*a)*d1^2*d2^2+(2*a^2+4*a+4)*d1*d2^3
  //        +(4*a^2+1)*d2^4+(2*a^2+3*a+1)*d1^3+(2*a^2+2*a+3)*d1^2*d2+(4*a^2+3*a+1)*d1*d2^2
  //        +(4*a^2+2*a)*d2^3+(4*a+1)*d1^2+(4*a^2+3)*d1*d2+(3*a^2+2*a+3)*d2^2
  int tempCoeffsf1f2[25*DEGREEEXTENSION] = {0,0,0, 0,0,0, 1,4,0, 1,3,2, 0,3,2,
                                            0,0,0, 3,0,4, 3,2,2, 3,4,4, 0,0,0,
                                            3,2,3, 1,3,4, 0,4,3, 0,0,0, 0,0,0,
                                            0,2,4, 4,4,2, 0,0,0, 0,0,0, 0,0,0,
                                            1,0,4, 0,0,0, 0,0,0, 0,0,0, 0,0,0};
  
  struct OrePoly *f1f2 =
    getOrePolyViaIntegerCoefficients(4,4,&Hom1,&Hom2, tempCoeffsf1f2);
  //1.a
  if (!isZero_OrePoly(z))
  {
    printf("1.a failed. Zero polynomial not recognized as being\
 zero\n");
    testSuccess = 0;
  }
  //1.b
  if (isZero_OrePoly(f1))
  {
    printf("1.b failed. f1 is recognized of being zero.\n");
    testSuccess = 0;
  }
  //2.a
  if (!isEqual_OrePoly(f2,f2))
  {
    printf("2.a failed. Two equal polynomials not recognized as\
 such.\n");
    testSuccess = 0;
  }
  //2.b
  if (isEqual_OrePoly(f1,f2))
  {
    printf("2.b failed. Two non-equal polynomials recognized as\
 equal.\n");
    testSuccess = 0;
  }
  //3.a
  struct OrePoly *k = scalarMult(0,f1);
  if (!isZero_OrePoly(k))
  {
    testSuccess = 0;
    tempOutp = OrePolyToString(k);
    printf("3.a failed. Scalarmultiplication with 0 did not result in\
 zero polynomial, but in: %s.\n", tempOutp);
    free(tempOutp);
  }
  free(k->coeffs);
  free(k);

  //3.b
  k = scalarMult(1,f2);
  if (!isEqual_OrePoly(k,f2))
  {
    testSuccess = 0;
    tempOutp = OrePolyToString(k);
    printf("3.b failed. Scalarmultiplication with 1 did not result in\
 the same polynomial f2, but in: %s.\n", tempOutp);
    free(tempOutp);
  }
  free(k->coeffs);
  free(k);

  //3.c
  k = scalarMult(2,f1);
  int i;
  for (i = 0; i<(k->degD1+1)*(k->degD2+1); ++i)
  {
    if (!isEqual_GF(k->coeffs[i], scalarMultGF(2,f1->coeffs[i])))
    {
      testSuccess = 0;
      tempOutp = OrePolyToString(k);
      printf("3.c failed. Scalarmultiplication with 2 did not result in\
 the polynomial 2*f1, but in: %s.\n", tempOutp);
      free(tempOutp);
      break;
    }
  }
  free(k->coeffs);
  free(k);

  //4.a
  k = add(z,z);
  if (!isZero_OrePoly(k))
  {
    testSuccess = 0;
    tempOutp = OrePolyToString(k);
    printf("4.a failed. Adding zero to zero results in: %s.\n", tempOutp);
    free(tempOutp);
  }
  free(k->coeffs);
  free(k);

  //4.b
  k = add(z, f3);
  if (!isEqual_OrePoly(k, f3))
  {
    testSuccess = 0;
    tempOutp = OrePolyToString(k);
    printf("4.b failed. Adding zero to f3 results in: %s.\n", tempOutp);
    free(tempOutp);
  }
  free(k->coeffs);
  free(k);

  //4.c
  k = add (f2,f3);
  for (i = 0; i<(k->degD1+1)*(k->degD2+1); ++i)
  {
    if (!isEqual_GF(k->coeffs[i],addGF(f2->coeffs[i],f3->coeffs[i])))
    {
      testSuccess = 0;
      tempOutp = OrePolyToString(k);
      printf("4.c failed. f2+f3 results in: %s.\n", tempOutp);
      free(tempOutp);
      break;
    }
  }
  free(k->coeffs);
  free(k);
  
  //5.a
  k = mult(z,f1);
  if(!isZero_OrePoly(k))
  {
    testSuccess = 0;
    tempOutp = OrePolyToString(k);
    printf("5.a failed. 0*f1 results in: %s.\n", tempOutp);
    free(tempOutp);
  }
  free(k->coeffs);
  free(k);

  //5.b
  k = mult(z,z);
  if(!isZero_OrePoly(k))
  {
    testSuccess = 0;
    tempOutp = OrePolyToString(k);
    printf("5.a failed. 0*0 results in: %s.\n", tempOutp);
    free(tempOutp);
  }
  free(k->coeffs);
  free(k);

  //5.c
  k = mult(one,one);
  if(!isEqual_OrePoly(k,one))
  {
    testSuccess = 0;
    tempOutp = OrePolyToString(k);
    printf("5.c failed. 1*1 results in: %s.\n", tempOutp);
    free(tempOutp);
  }
  free(k->coeffs);
  free(k);

  //5.d
  k = mult(one,f1);
  if(!isEqual_OrePoly(k,f1))
  {
    testSuccess = 0;
    tempOutp = OrePolyToString(k);
    printf("5.d failed. 1*f1 results in: %s.\n", tempOutp);
    free(tempOutp);
  }
  free(k->coeffs);
  free(k);

  //5.e
  k = mult(f1,f2);
  if(!isEqual_OrePoly(k,f1f2))
  {
    testSuccess = 0;
    tempOutp = OrePolyToString(k);
    printf("5.e failed. f1*f2 results in: %s.\n", tempOutp);
    free(tempOutp);
  }
  free(k->coeffs);
  free(k);

  //freeing memory and returning
  free(z->coeffs);
  free(z);
  free(one->coeffs);
  free(one);
  free(f1->coeffs);
  free(f1);
  free(f2->coeffs);
  free(f2);
  free(f3->coeffs);
  free(f3);
  free(f1f2->coeffs);
  free(f1f2);
  return testSuccess;
}//test_OrePoly_Arith

/**
   This tests if the polynomials that are randomly generated
   and which are supposed to commute are in fact commuting.
   The following cases are covered:
   1. Randomly generate polynomials from degree 0 - 5 and check for commutation with itself
   2. Use the generated polynomials from 1. and create keys from those and check for commutation (degree 2-6).
 */
int test_SecretKey_Validity()
{
  int testSuccess = 1;
  char *tempOutp;

  struct OrePoly *f1 = getRandomOrePoly(0,0,&Hom1,&Hom2);
  struct OrePoly *f2 = getRandomOrePoly(1,0,&Hom1,&Hom2);
  struct OrePoly *f3 = getRandomOrePoly(1,1,&Hom1,&Hom2);
  struct OrePoly *f4 = getRandomOrePoly(1,2,&Hom1,&Hom2);
  struct OrePoly *f5 = getRandomOrePoly(2,2,&Hom1,&Hom2);
  struct OrePoly *f6 = getRandomOrePoly(3,2,&Hom1,&Hom2);

  //1.
  struct OrePoly *k1 = mult(f1,f1);
  struct OrePoly *k2 = minus(k1,k1);
  if(!isZero_OrePoly(k2))
  {
    testSuccess = 0;
    tempOutp = OrePolyToString(k2);
    printf("Test 1.0 failed. f1*f1 - f1*f1 is not equal to 0, but to: %s\n", tempOutp);
    free(tempOutp);
  }
  free(k1->coeffs);
  free(k1);
  free(k2->coeffs);
  free(k2);

  k1 = mult (f2,f2);
  k2 = minus(k1,k1);
  if(!isZero_OrePoly(k2))
  {
    testSuccess = 0;
    tempOutp = OrePolyToString(k2);
    printf("Test 1.1 failed. f2*f2 - f2*f2 is not equal to 0, but to: %s\n", tempOutp);
    free(tempOutp);
  }
  free(k1->coeffs);
  free(k1);
  free(k2->coeffs);
  free(k2);

  k1 = mult (f3,f3);
  k2 = minus(k1,k1);
  if(!isZero_OrePoly(k2))
  {
    testSuccess = 0;
    tempOutp = OrePolyToString(k2);
    printf("Test 1.2 failed. f3*f3 - f3*f3 is not equal to 0, but to: %s\n", tempOutp);
    free(tempOutp);
  }
  free(k1->coeffs);
  free(k1);
  free(k2->coeffs);
  free(k2);

  k1 = mult (f4,f4);
  k2 = minus(k1,k1);
  if(!isZero_OrePoly(k2))
  {
    testSuccess = 0;
    tempOutp = OrePolyToString(k2);
    printf("Test 1.3 failed. f4*f4 - f4*f4 is not equal to 0, but to: %s\n", tempOutp);
    free(tempOutp);
  }
  free(k1->coeffs);
  free(k1);
  free(k2->coeffs);
  free(k2);

  k1 = mult (f5,f5);
  k2 = minus(k1,k1);
  if(!isZero_OrePoly(k2))
  {
    testSuccess = 0;
    tempOutp = OrePolyToString(k2);
    printf("Test 1.4 failed. f5*f5 - f5*f5 is not equal to 0, but to: %s\n", tempOutp);
    free(tempOutp);
  }
  free(k1->coeffs);
  free(k1);
  free(k2->coeffs);
  free(k2);

  k1 = mult (f6,f6);
  k2 = minus(k1,k1);
  if(!isZero_OrePoly(k2))
  {
    testSuccess = 0;
    tempOutp = OrePolyToString(k2);
    printf("Test 1.5 failed. f6*f6 - f6*f6 is not equal to 0, but to: %s\n", tempOutp);
    free(tempOutp);
  }
  free(k1->coeffs);
  free(k1);
  free(k2->coeffs);
  free(k2);

  //2.
  k1 = generateRandomSecretKey(2,f1);
  k2 = generateRandomSecretKey(3,f1);
  struct OrePoly *k3 = mult(k1,k2);
  struct OrePoly *k4 = mult(k2,k1);
  struct OrePoly *k5 = minus(k3,k4);
  if (!isZero_OrePoly(k5))
  {
    testSuccess = 0;
    tempOutp = OrePolyToString(k5);
    printf("Test 2.1 failed. Generated random secret keys did not commute. Their difference is: %s\n", tempOutp);
    free(tempOutp);
    tempOutp = OrePolyToString(k1);
    printf("The first polynomial was: %s\n", tempOutp);
    free(tempOutp);
    tempOutp = OrePolyToString(k2);
    printf("The second polynomial was: %s\n", tempOutp);
    free(tempOutp);
  }
  free(k1->coeffs);
  free(k1);
  free(k2->coeffs);
  free(k2);
  free(k3->coeffs);
  free(k3);
  free(k4->coeffs);
  free(k4);
  free(k5->coeffs);
  free(k5);

  k1 = generateRandomSecretKey(5,f2);
  k2 = generateRandomSecretKey(4,f2);
  k3 = mult(k1,k2);
  k4 = mult(k2,k1);
  k5 = minus(k3,k4);
  if (!isZero_OrePoly(k5))
  {
    testSuccess = 0;
    tempOutp = OrePolyToString(k5);
    printf("Test 2.2 failed. Generated random secret keys did not commute. Their difference is: %s\n", tempOutp);
    free(tempOutp);
    tempOutp = OrePolyToString(k1);
    printf("The first polynomial was: %s\n", tempOutp);
    free(tempOutp);
    tempOutp = OrePolyToString(k2);
    printf("The second polynomial was: %s\n", tempOutp);
    free(tempOutp);
  }
  free(k1->coeffs);
  free(k1);
  free(k2->coeffs);
  free(k2);
  free(k3->coeffs);
  free(k3);
  free(k4->coeffs);
  free(k4);
  free(k5->coeffs);
  free(k5);
  
  k1 = generateRandomSecretKey(3,f3);
  k2 = generateRandomSecretKey(6,f3);
  k3 = mult(k1,k2);
  k4 = mult(k2,k1);
  k5 = minus(k3,k4);
  if (!isZero_OrePoly(k5))
  {
    testSuccess = 0;
    tempOutp = OrePolyToString(k5);
    printf("Test 2.3 failed. Generated random secret keys did not commute. Their difference is: %s\n", tempOutp);
    free(tempOutp);
    tempOutp = OrePolyToString(k1);
    printf("The first polynomial was: %s\n", tempOutp);
    free(tempOutp);
    tempOutp = OrePolyToString(k2);
    printf("The second polynomial was: %s\n", tempOutp);
    free(tempOutp);
  }
  free(k1->coeffs);
  free(k1);
  free(k2->coeffs);
  free(k2);
  free(k3->coeffs);
  free(k3);
  free(k4->coeffs);
  free(k4);
  free(k5->coeffs);
  free(k5);

  k1 = generateRandomSecretKey(5,f4);
  k2 = generateRandomSecretKey(6,f4);
  k3 = mult(k1,k2);
  k4 = mult(k2,k1);
  k5 = minus(k3,k4);
  if (!isZero_OrePoly(k5))
  {
    testSuccess = 0;
    tempOutp = OrePolyToString(k5);
    printf("Test 2.4 failed. Generated random secret keys did not commute. Their difference is: %s\n", tempOutp);
    free(tempOutp);
    tempOutp = OrePolyToString(k1);
    printf("The first polynomial was: %s\n", tempOutp);
    free(tempOutp);
    tempOutp = OrePolyToString(k2);
    printf("The second polynomial was: %s\n", tempOutp);
    free(tempOutp);
  }
  free(k1->coeffs);
  free(k1);
  free(k2->coeffs);
  free(k2);
  free(k3->coeffs);
  free(k3);
  free(k4->coeffs);
  free(k4);
  free(k5->coeffs);
  free(k5);

  k1 = generateRandomSecretKey(6,f5);
  k2 = generateRandomSecretKey(6,f5);
  k3 = mult(k1,k2);
  k4 = mult(k2,k1);
  k5 = minus(k3,k4);
  if (!isZero_OrePoly(k5))
  {
    testSuccess = 0;
    tempOutp = OrePolyToString(k5);
    printf("Test 2.5 failed. Generated random secret keys did not commute. Their difference is: %s\n", tempOutp);
    free(tempOutp);
    tempOutp = OrePolyToString(k1);
    printf("The first polynomial was: %s\n", tempOutp);
    free(tempOutp);
    tempOutp = OrePolyToString(k2);
    printf("The second polynomial was: %s\n", tempOutp);
    free(tempOutp);
  }
  free(k1->coeffs);
  free(k1);
  free(k2->coeffs);
  free(k2);
  free(k3->coeffs);
  free(k3);
  free(k4->coeffs);
  free(k4);
  free(k5->coeffs);
  free(k5);

  k1 = generateRandomSecretKey(3,f6);
  k2 = generateRandomSecretKey(3,f6);
  k3 = mult(k1,k2);
  k4 = mult(k2,k1);
  k5 = minus(k3,k4);
  if (!isZero_OrePoly(k5))
  {
    testSuccess = 0;
    tempOutp = OrePolyToString(k5);
    printf("Test 2.6 failed. Generated random secret keys did not commute. Their difference is: %s\n", tempOutp);
    free(tempOutp);
    tempOutp = OrePolyToString(k1);
    printf("The first polynomial was: %s\n", tempOutp);
    free(tempOutp);
    tempOutp = OrePolyToString(k2);
    printf("The second polynomial was: %s\n", tempOutp);
    free(tempOutp);
  }
  free(k1->coeffs);
  free(k1);
  free(k2->coeffs);
  free(k2);
  free(k3->coeffs);
  free(k3);
  free(k4->coeffs);
  free(k4);
  free(k5->coeffs);
  free(k5);

  free(f1->coeffs);
  free(f1);
  free(f2->coeffs);
  free(f2);
  free(f3->coeffs);
  free(f3);
  free(f4->coeffs);
  free(f4);
  free(f5->coeffs);
  free(f5);
  free(f6->coeffs);
  free(f6);
  return testSuccess;
}

/**
   This is a testsuite to check if the elements in the Galois field are
   correctly enumerated by the function getAllPossibleGFElements.
   The tests cover random checks, namely:
   1. The first element is 0
   2. The last element is 4a^(DEGREEEXTENSION-1) + ... + 4a+ 4
   3. The element at position 38 ist 1 + 2a + 3a^2 %TODO: Alter this test to work in a general framework
 */
int test_GF_Enumeration()
{//test_GF_Enumeration
  struct GFModulus *obj = getAllPossibleGFElements();
  int testSuccess = 1;
  char* tempOutp;
  if (obj[0].coeffs[0] != 0 || obj[0].coeffs[1] != 0 || obj[0].coeffs[2] !=0)
  {
    tempOutp = GFModulusToString(obj[0]);
    printf("Test 1 failed. Expected the zero-element at the 0-position, but got: %s\n", tempOutp);
    free(tempOutp);
    testSuccess = 0;
  }
  if (obj[124].coeffs[0] != 4 || obj[124].coeffs[1] != 4 || obj[124].coeffs[2] !=4)
  {
    tempOutp = GFModulusToString(obj[124]);
    printf("Test 2 failed. Expected the element 4 + 4a +... at the 124-position, but got: %s\n", tempOutp);
    free(tempOutp);
    testSuccess = 0;
  }
  if (obj[38].coeffs[0] != 1 || obj[38].coeffs[1] != 2 || obj[38].coeffs[2] !=3)
  {
    tempOutp = GFModulusToString(obj[38]);
    printf("Test 3 failed. Expected the element 1 + 2a +... at the 38-position, but got: %s\n", tempOutp);
    free(tempOutp);
    testSuccess = 0;
  }

  /* int i; */
  /* for (i = 0; i<NUMBEROFELEMENTSINGF; ++i) */
  /* { */
  /*   tempOutp = GFModulusToString(obj[i]); */
  /*   printf("%s\n", tempOutp); */
  /*   free(tempOutp); */
  /* } */
  free(obj);
  return testSuccess;
}//test_GF_Enumeration

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
  if (test_OrePoly_Arith() == 0)
    printf("test_OrePoly_Arith failed.\n");
  else
    printf("test_OrePoly_Arith succeeded.\n");
  if (test_SecretKey_Validity() == 0)
    printf("test_SecretKey_Validity failed.\n");
  else
    printf("test_SecretKey_Validity succeeded.\n");
  if (test_GF_Enumeration() ==0)
    printf("test_GF_Enumeration failed.\n");
  else
    printf("test_GF_Enumeration succeeded.\n");
  return 0;
}//main
