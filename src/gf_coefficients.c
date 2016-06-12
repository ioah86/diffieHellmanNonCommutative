#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "gf_coefficients.h"


//////////////////////////////////////////////////
// CAUTION: The Following functions and definitions have to be made
// mathematically correct; i.e. a person who modifies this part of the
// code has to know what he/she is doing.
// We assume the minimal polynomial of a to be
// a^3 + 3*a + 3
//////////////////////////////////////////////////

char* GFModulusToString(struct GFModulus inp)
{
  int coeffInStringSize = MODULUS/10+1;
  int exponentInStringSize = DEGREEEXTENSION/10+1;
  char* result =
    malloc(DEGREEEXTENSION*(coeffInStringSize+5+exponentInStringSize));
  int i;
  char tmpStr[coeffInStringSize+6+exponentInStringSize];
  strcpy(result, "(");
  for (i = 0; i<DEGREEEXTENSION; ++i)
  {
    if (i!=DEGREEEXTENSION-1)
      sprintf(tmpStr,"%ia^%i + ",inp.coeffs[i],i);
    else
      sprintf(tmpStr,"%ia^%i",inp.coeffs[i],i);
    strcat(result,tmpStr);
  }
  strcat(result, ")");
  return result;
}

void GFModulusToStdOut(struct GFModulus inp)
{//GFModulusToStdOut
  int i;
  printf("(");
  for (i=0; i<DEGREEEXTENSION;++i)
  {
    if (i!=DEGREEEXTENSION-1)
      printf("%ia^%i + ",inp.coeffs[i],i);
    else
      printf("%ia^%i)",inp.coeffs[i],i);
  }
}//GFModulusToStdOut

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
    int a_0 = inp1.coeffs[0];
  int a_1 = inp1.coeffs[1];
  int a_2 = inp1.coeffs[2];
  int b_0 = inp2.coeffs[0];
  int b_1 = inp2.coeffs[1];
  int b_2 = inp2.coeffs[2];
  result.coeffs[0] = (a_0*b_0 + 2*a_2*b_1 + 2*a_1*b_2) %MODULUS;
  result.coeffs[1] = (a_1*b_0 + a_0*b_1 + 2*a_2*b_1 + 2*a_1*b_2 + 2*a_2*b_2) %MODULUS;
  result.coeffs[2] = (a_2*b_0 + a_1*b_1 + a_0*b_2 + 2*a_2*b_2) %MODULUS;

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

//Homomorphisms begin
struct GFModulus Hom1(struct GFModulus inp)
{
  //Short description of Hom1:
  //a |--> 3*a^2 + 1
  struct GFModulus result;
  int a_0 = inp.coeffs[0];
  int a_1 = inp.coeffs[1];
  int a_2 = inp.coeffs[2];
  result.coeffs[0] = (a_0 + a_1 + a_2) %MODULUS;
  result.coeffs[1] = (3*a_2) %MODULUS;
  result.coeffs[2] = (3*a_1 + 4*a_2) %MODULUS;
  return result;
}

struct GFModulus Hom2(struct GFModulus inp)
{
  //Short description of Hom2:
  //a |--> 2*a^2 + 4*a + 4
  struct GFModulus result;
  int a_0 = inp.coeffs[0];
  int a_1 = inp.coeffs[1];
  int a_2 = inp.coeffs[2];
  result.coeffs[0] = (a_0 + 4*a_1 + 3*a_2) %MODULUS;
  result.coeffs[1] = (4*a_1 + 2*a_2) %MODULUS;
  result.coeffs[2] = (2*a_1) %MODULUS;
  return result;
}


//Homomorphisms END

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
{//getAllPossibleGFElements
  struct GFModulus *result = malloc(NUMBEROFELEMENTSINGF*sizeof(struct GFModulus));
  long i; long tempi;
  int j;
  for (i = 0; i<NUMBEROFELEMENTSINGF; ++i)
  {//iterate through all elements in the GF
    tempi = i;
    for (j = 0; j<DEGREEEXTENSION; ++j)
    {//iterate through all coefficients
        result[i].coeffs[j] = tempi % MODULUS;
	tempi/=MODULUS;
    }//iterate through all coefficients
  }//iterate through all elements in the GF
  return result;
}//getAllPossibleGFElements
