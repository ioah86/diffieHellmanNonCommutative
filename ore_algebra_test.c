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
 5. testing the multiplication TODO
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
  //f1 =(a^2+a+3)*d1^2+(4*a^2+4)*d1*d2+(a^2+4*a)*d2^2+(3*a^2+4a+4)*d1+(a^2+3)*d2
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
  return testSuccess;
}//test_OrePoly_Arith

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
  return 0;
}//main
