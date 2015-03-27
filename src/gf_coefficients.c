#include <stdlib.h>
#include <stdio.h>
#include "gf_coefficients.h"


//////////////////////////////////////////////////
// CAUTION: The Following functions and definitions have to be made
// mathematically correct; i.e. a person who modifies this part of the
// code has to know what he/she is doing.
// We assume the minimal polynomial of a to be
// x^3 + 3*x + 3
//////////////////////////////////////////////////

char* GFModulusToString(struct GFModulus inp)
{
  int coeffInStringSize = MODULUS/10+1;
  char* result = malloc(DEGREEEXTENSION*coeffInStringSize + 13);
  sprintf(result,"(%i + %ia + %ia^2)",inp.coeffs[0],inp.coeffs[1],
		  inp.coeffs[2]);
  return result;
}

void GFModulusToStdOut(struct GFModulus inp)
{
  printf("(%i + %ia + %ia^2)",inp.coeffs[0],inp.coeffs[1],inp.coeffs[2]);
}

struct GFModulus addGF(struct GFModulus inp1,struct GFModulus inp2)
{//addGF
  struct GFModulus result;
  int i;
  for (i = 0; i<=DEGREEEXTENSION; ++i)
  {
    result.coeffs[i] = (inp1.coeffs[i] + inp2.coeffs[i]) % MODULUS;
  }
  return result;
}//addGF

struct GFModulus minusGF(struct GFModulus inp1, struct GFModulus inp2)
{//minusGF
  struct GFModulus result;
  int i;
  for (i = 0; i<=DEGREEEXTENSION; ++i)
  {
    result.coeffs[i] = (inp1.coeffs[i] - inp2.coeffs[i]) % MODULUS;
    if (result.coeffs[i]<0)
      result.coeffs[i] += MODULUS;
  }
  return result;
}//minusGF

struct GFModulus multGF(struct GFModulus inp1, struct GFModulus inp2)
{//multGF
  struct GFModulus result;
  int b,c,d,e,f,g;
  b = inp1.coeffs[0];
  c = inp1.coeffs[1];
  d = inp1.coeffs[2];
  e = inp2.coeffs[0];
  f = inp2.coeffs[1];
  g = inp2.coeffs[2];
  result.coeffs[0] = (b*e + 2*d*f + 2*c*g) %MODULUS;
  result.coeffs[1] = (c*e + b*f + 2*d*f + 2*c*g + 2*d*g)%MODULUS;
  result.coeffs[2] = (d*e + c*f + b*g + 2*d*g)%MODULUS;
  return result;
}//multGF

int isZero_GF(struct GFModulus inp)
{//isZero_GF
  int i;
  for (i = 0; i<DEGREEEXTENSION; ++i)
  {
    if (inp.coeffs[i] %MODULUS !=0)
    {
      return 0;
    }
  }
  return 1;
}//isZero_GF

int isEqual_GF(struct GFModulus inp1, struct GFModulus inp2)
{//isEqual_OrePoly
  int i;
  for (i=0;i<DEGREEEXTENSION; ++i)
  {
    if (inp1.coeffs[i] %MODULUS != inp2.coeffs[i] %MODULUS)
    {
      return 0;
    }
  }
  return 1;
}//isEqual_OrePoly

struct GFModulus getZeroElemGF()
{//getZeroElemGF
  struct GFModulus result;
  result.coeffs[0]=0;
  int i;
  for (i = 0; i<DEGREEEXTENSION; ++i)
    result.coeffs[i] = 0;
  return result;
}//getZeroElemGF

struct GFModulus getIdentityElemGF()
{//getIdentityElemGF
  struct GFModulus result;
  result.coeffs[0] = 1;
  int i;
  for (i = 1; i<DEGREEEXTENSION; ++i)
    result.coeffs[i] = 0;
  return result;
}//getIdentityElemGF


struct GFModulus getMinusOneElemGF()
{//getMinusOneElemGF
  struct GFModulus result = getIdentityElemGF();
  result.coeffs[0] = MODULUS -1;
  return result;
}//getMinusOneElemGF

struct GFModulus getRandomGFElem()
{//getRandomGFElem
  struct GFModulus result;
  int i;
  for (i = 0; i<DEGREEEXTENSION; ++i)
    result.coeffs[i] = rand() % MODULUS;
  return result;
}//getRandomGFElem


struct GFModulus identityMap(struct GFModulus inp)
{//identityMap
  return inp;
}//identityMap

struct GFModulus Hom1(struct GFModulus inp)
{//Hom1
  //Short description of Hom1:
  //a |--> 3*a^2 + 1
  int b,c,d;
  b = inp.coeffs[0];
  c = inp.coeffs[1];
  d = inp.coeffs[2];
  struct GFModulus result;
  result.coeffs[0] = (b + c + d)%MODULUS;
  result.coeffs[1] = (3*d) %MODULUS;
  result.coeffs[2] = (3*c + 4*d) %MODULUS;
  return result;
}//Hom1

struct GFModulus Hom2(struct GFModulus inp)
{//Hom2
  //Short description of Hom2:
  //a |--> 2*a^2 + 4*a + 4
  int b,c,d;
  b = inp.coeffs[0];
  c = inp.coeffs[1];
  d = inp.coeffs[2];
  struct GFModulus result;
  result.coeffs[0] = (b + 4*c + 3*d)%MODULUS;
  result.coeffs[1] = (4*c + 2*d) %MODULUS;
  result.coeffs[2] = (2*c) %MODULUS;
  return result;
}//Hom2


struct GFModulus scalarMultGF(int s, struct GFModulus inp)
{//scalarMultGF
  struct GFModulus result;
  int scalarNormed = s %MODULUS;
  if (scalarNormed < 0)
    scalarNormed += MODULUS;
  int i;
  for (i = 0; i<DEGREEEXTENSION; ++i)
    result.coeffs[i] = (scalarNormed*inp.coeffs[i])%MODULUS;
  return result;
}//scalarMultGF

//////////////////////////////////////////////////
//CAUTION FINISHED
//////////////////////////////////////////////////
struct GFModulus* getAllPossibleGFElements()
{//getAllPossibleGFElements()
  //TODO: Make that work for arbitrary fields
  struct GFModulus *result = malloc(NUMBEROFELEMENTSINGF*sizeof(struct GFModulus));
  int i; int j; int k;
  int l = 0;
  for (i = 0; i<MODULUS; ++i)
  {
    for (j = 0; j<MODULUS; ++j)
    {
      for (k = 0; k<MODULUS; ++k)
      {
        result[l].coeffs[0] = i;
        result[l].coeffs[1] = j;
        result[l].coeffs[2] = k;
        ++l;
      }
    }
  }
  return result;
}//getAllPossibleGFElements
