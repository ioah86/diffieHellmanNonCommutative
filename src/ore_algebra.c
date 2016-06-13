/*
  Implementation of the Diffie-Hellman protocol -- Non-commutative version

  (c) Reinhold Burger and Albert Heinle
 */

#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#include "ore_algebra.h"


char* OrePolyToString(struct OrePoly* inp)
{//OrePolyToString
  int i; int j;//iteration variables
  int coeffInStringSize = DEGREEEXTENSION*(MODULUS/10+7+DEGREEEXTENSION/10);
  int monomialInStringSize = (inp->degD1)/10+(inp->degD2)/10 + 8;
  char* coeffInString; // encoding of the coefficient as
				     // string
  char monomialInString[monomialInStringSize];
  int alreadyACoeffFound = 0;
  int resultSize=
    (inp->degD1+1)*(inp->degD2+1)*
    (coeffInStringSize + monomialInStringSize +3);
  char* result = malloc(resultSize * sizeof(char));
  if (result ==NULL)
  {
    printf("Error allocating memory for the string representation of OrePoly.\n");
    exit(1);
  }
  result[0] = '\0';
  for(i=0; i<=inp->degD2; ++i)
  {//iterating through all possible degrees of D2
    for (j = 0; j<= inp->degD1; ++j)
    {//iterating through all possible degrees of D1
      if (isZero_GF(inp->coeffs[i*(inp->degD1+1)+j])==0)
      {//only then, we are printing something
	if (alreadyACoeffFound)
	  strcat(result, " + ");
	else
	  alreadyACoeffFound = 1;
	coeffInString = GFModulusToString((inp->coeffs)[i*(inp->degD1+1) + j]);
        strcat(result,coeffInString);
	free(coeffInString);
	if (i !=0 || j!=0)
	{
	  if (i !=0 && j!=0)
	  {
	    sprintf(monomialInString,"d1^%id2^%i", j, i);
	    strcat(result, monomialInString);
	  }
	  else if (i==0)
	  {
	    sprintf(monomialInString,"d1^%i", j);
	    strcat(result, monomialInString);
	  }
	  else
	  {
	    sprintf(monomialInString,"d2^%i", i);
	    strcat(result, monomialInString);
	  }
	}
      }//only then, we are printing something
    }//iterating through all possible degrees of D1
  }//iterating through all possible degrees of D2
  if(!alreadyACoeffFound)
  {
    result[0] = '0';
    result[1] = '\0';
  }
  return result;
}//OrePolyToString

void OrePolyToStdOut(struct OrePoly* inp)
{//OrePolyToStdOut
  int i; int j;//iteration variables
  int alreadyACoeffFound = 0;
  for(i=0; i<=inp->degD2; ++i)
  {//iterating through all possible degrees of D2
    for (j = 0; j<= inp->degD1; ++j)
    {//iterating through all possible degrees of D1
      if (isZero_GF(inp->coeffs[i*(inp->degD1+1)+j])==0)
      {//only then, we are printing something
	if (alreadyACoeffFound)
	  printf(" + ");
	else
	  alreadyACoeffFound = 1;
	GFModulusToStdOut((inp->coeffs)[i*(inp->degD1+1) + j]);
        if (i !=0 || j!=0)
	{
	  if (i !=0 && j!=0)
	  {
	    printf("d1^%id2^%i", j, i);
	  }
	  else if (i==0)
	  {
	    printf("d1^%i", j);
	  }
	  else
	  {
	    printf("d2^%i", i);
	  }
	}
      }//only then, we are printing something
    }//iterating through all possible degrees of D1
  }//iterating through all possible degrees of D2
  if(!alreadyACoeffFound)
  {
    printf("0");
  }
  printf("\n");
}//OrePolyToStdOut

struct OrePoly * getOrePolyViaIntegerCoefficients(int degD1, int degD2,
						struct GFModulus
						(*ptrD1manip)(struct GFModulus),
						struct	GFModulus
						(*ptrD2manip)(struct GFModulus),
						int *coeffs)
{//getOrePolyViaIntegerCoefficients
  struct OrePoly *result = malloc(sizeof(struct OrePoly));
  result->degD1 = degD1;
  result->degD2 = degD2;
  result->ptrD1manip = ptrD1manip;
  result->ptrD2manip = ptrD2manip;
  int i;
  int j = 0;
  int k;
  struct GFModulus *resCoeffs = malloc(sizeof(struct GFModulus)*(degD1+1)*(degD2+1));
  for (i = 0; i<(degD1+1)*(degD2+1)*DEGREEEXTENSION;
       i+=DEGREEEXTENSION)
  {
    for (k = 0; k<DEGREEEXTENSION; ++k)
    {
      resCoeffs[j].coeffs[k] = coeffs[i+k] % MODULUS;
      if (resCoeffs[j].coeffs[k] < 0)
	resCoeffs[j].coeffs[k] += MODULUS;
    }
    j+=1;
  }
  result->coeffs = resCoeffs;
  return result;
}//getOrePolyViaIntegerCoefficients

int isZero_OrePoly(struct OrePoly *inp)
{//isZero_OrePoly
  int i;
  for (i=0; i<(inp->degD1+1)*(inp->degD2+1); ++i)
  {
    if (isZero_GF(inp->coeffs[i])==0)
    {
      return 0;
    }
  }
  return 1;
}//isZero_OrePoly

int isEqual_OrePoly(struct OrePoly *inp1, struct OrePoly *inp2)
{//isEqual_OrePoly
  if (inp1->degD1 != inp2->degD1)
    return 0;
  if (inp1->degD2 != inp2->degD2)
    return 0;
  if (inp1->ptrD1manip != inp2->ptrD1manip)
    return 0;
  if (inp1->ptrD2manip != inp2->ptrD2manip)
    return 0;
  int i; int j;
  for (i = 0; i<(inp1->degD1+1)*(inp1->degD2+1); ++i)
  {
    for (j=0; j<DEGREEEXTENSION; ++j)
    {
      if (inp1->coeffs[i].coeffs[j] != inp2->coeffs[i].coeffs[j])
      {
        return 0;
      }
    }
  }
  return 1;
}//isEqual_OrePoly

struct OrePoly *getIdentityElemOrePoly(struct GFModulus
				       (*ptrD1manip)(struct GFModulus),
				       struct GFModulus
				       (*ptrD2manip)(struct GFModulus))
{//getIdentityElemOrePoly
  struct OrePoly *result = malloc(sizeof(struct OrePoly));
  result->degD1 = 0;
  result->degD2 = 0;
  result->ptrD1manip = ptrD1manip;
  result->ptrD2manip = ptrD2manip;
  struct GFModulus *coeffs = malloc(sizeof(struct GFModulus));
  coeffs[0] = getIdentityElemGF();
  result->coeffs = coeffs;
  return result;
}//getIdentityElemOrePoly

struct OrePoly *getZeroElemOrePoly(struct GFModulus
				   (*ptrD1manip)(struct GFModulus),
				   struct GFModulus
				   (*ptrD2manip)(struct GFModulus))
{//getZeroElemOrePoly
  struct OrePoly *result = malloc(sizeof(struct OrePoly));
  result->degD1 = 0;
  result->degD2 = 0;
  result->ptrD1manip = ptrD1manip;
  result->ptrD2manip = ptrD2manip;
  struct GFModulus *coeffs = malloc(sizeof(struct GFModulus));
  coeffs[0] = getZeroElemGF();
  result->coeffs = coeffs;
  return result;
}//getZeroElemOrePoly

struct OrePoly* scalarMult(int s, struct OrePoly* inp)
{//scalarMultipy
  int i;
  struct OrePoly* result = malloc(sizeof(struct OrePoly));
  if (result==NULL)
  {
    printf("Could not allocate memory for the result of the \
scalar multiplication.\n");
    exit(1);
  }
  if (s % MODULUS == 0)
  {//multiplying by 0
    result->degD1 = 0;
    result->degD2 = 0;
    struct GFModulus * coeffs = malloc(sizeof(struct GFModulus));
    if (coeffs == NULL)
    {
      printf("Could not reserve memory for one integer in \
scalarMult.\n");
      exit(1);
    }
    struct GFModulus zeroElem = getZeroElemGF();
    coeffs[0] = zeroElem;
    result->coeffs = coeffs;
    result->ptrD1manip = inp->ptrD1manip;
    result->ptrD2manip = inp->ptrD2manip;
    return result;
  }//multiplying by 0
  result->degD1 = inp->degD1;
  result->degD2 = inp->degD2;
  struct GFModulus *coeffs =
    malloc((result->degD1+1)*(result->degD2+1)*sizeof(struct GFModulus));
  if (coeffs == NULL)
  {
    printf("Could not reserve memory for new OrePoly in \
ScalarMultiplication.\n");
    exit(1);
  }
  result->coeffs = coeffs;
  for (i = 0; i<(result->degD1+1)*(result->degD2+1); ++i)
  {
    result->coeffs[i] = (scalarMultGF(s,(inp->coeffs[i])));
  }
  result->ptrD1manip = inp->ptrD1manip;
  result->ptrD2manip = inp->ptrD2manip;
  return result;
}//scalarMultipy

struct OrePoly* add(struct OrePoly* inp1, struct OrePoly* inp2)
{//add
  int i; int j;
  if (inp1->ptrD1manip != inp2->ptrD1manip ||
      inp1->ptrD2manip != inp2->ptrD2manip)
  {//in this case, we try to add two different operators
    printf("Maps for the two operators we want to add do not coincide.\n");
    exit(1);
  }//in this case, we try to add two different operators
  struct OrePoly* result = malloc(sizeof(struct OrePoly));
  if (result==NULL)
  {
    printf("Could not allocate memory for the result of the \
addition.\n");
    exit(1);
  }
  if (inp1->degD1 > inp2->degD1)
    result->degD1 = inp1->degD1;
  else
    result->degD1 = inp2->degD1;
  if (inp1->degD2 > inp2->degD2)
    result->degD2 = inp1->degD2;
  else
    result->degD2 = inp2->degD2;
  struct GFModulus* coeffs =
    malloc((result->degD1+1)*(result->degD2+1)*sizeof(struct GFModulus));
  if (coeffs == NULL)
  {
    printf("Could not allocate memory for the coefficients of the \
result of the addition");
    exit(1);
  }
  //memset(coeffs, 0,
  //(result->degD1+1)*(result->degD2+1)*sizeof(int));
  for (i = 0; i<(result->degD1+1)*(result->degD2+1);++i)
    coeffs[i] = getZeroElemGF();
  result->coeffs = coeffs;
  result->ptrD1manip = inp1->ptrD1manip;
  result->ptrD2manip = inp1->ptrD2manip;
  for (i = 0; i <= inp1->degD2 ; ++i)
  {
    for (j=0; j<= inp1->degD1; ++j)
    {
      result->coeffs[i*(result->degD1+1) + j] =
	addGF(result->coeffs[i*(result->degD1+1) + j],
	      inp1->coeffs[i*(inp1->degD1+1) + j]);
    }
  }
  for (i = 0; i <= inp2->degD2 ; ++i)
  {
    for (j=0; j<= inp2->degD1; ++j)
    {
      result->coeffs[i*(result->degD1+1) + j] =
	addGF(result->coeffs[i*(result->degD1+1) + j],
	      inp2->coeffs[i*(inp2->degD1+1) + j]);
    }
  }
  return result;
}//add

struct OrePoly* minus(struct OrePoly* inp1, struct OrePoly* inp2)
{//minus
  int i; int j;
  if (inp1->ptrD1manip != inp2->ptrD1manip ||
      inp1->ptrD2manip != inp2->ptrD2manip)
  {//in this case, we try to minus two different operators
    printf("Maps for the two operators we want to minus do not coincide.\n");
    exit(1);
  }//in this case, we try to minus two different operators
  struct OrePoly* result = malloc(sizeof(struct OrePoly));
  if (result==NULL)
  {
    printf("Could not allocate memory for the result of the \
substraction.\n");
    exit(1);
  }
  if (inp1->degD1 > inp2->degD1)
    result->degD1 = inp1->degD1;
  else
    result->degD1 = inp2->degD1;
  if (inp1->degD2 > inp2->degD2)
    result->degD2 = inp1->degD2;
  else
    result->degD2 = inp2->degD2;
  struct GFModulus* coeffs =
    malloc((result->degD1+1)*(result->degD2+1)*sizeof(struct GFModulus));
  if (coeffs == NULL)
  {
    printf("Could not allocate memory for the coefficients of the \
result of the substraction");
    exit(1);
  }
  //memset(coeffs, 0,
  //(result->degD1+1)*(result->degD2+1)*sizeof(int));
  for (i = 0; i<(result->degD1+1)*(result->degD2+1);++i)
    coeffs[i] = getZeroElemGF();
  result->coeffs = coeffs;
  result->ptrD1manip = inp1->ptrD1manip;
  result->ptrD2manip = inp1->ptrD2manip;
  for (i = 0; i <= inp1->degD2 ; ++i)
  {
    for (j=0; j<= inp1->degD1; ++j)
    {
      result->coeffs[i*(result->degD1+1) + j] =
	addGF(result->coeffs[i*(result->degD1+1) + j],
	      inp1->coeffs[i*(inp1->degD1+1) + j]);
    }
  }
  for (i = 0; i <= inp2->degD2 ; ++i)
  {
    for (j=0; j<= inp2->degD1; ++j)
    {
      result->coeffs[i*(result->degD1+1) + j] =
	minusGF(result->coeffs[i*(result->degD1+1) + j],
	      inp2->coeffs[i*(inp2->degD1+1) + j]);
    }
  }
  return result;
}//minus

struct OrePoly* mult(struct OrePoly* inp1, struct OrePoly* inp2)
{//mult
  int i; int j; int k; int l; int m;
  struct GFModulus tempCoeff;
  if (inp1->ptrD1manip != inp2->ptrD1manip ||
      inp1->ptrD2manip != inp2->ptrD2manip)
  {//in this case, we try to multiply two different operators
    printf("Maps for the two operators we want to multiply do not coincide.\n");
    exit(1);
  }//in this case, we try to multiply two different operators
  struct OrePoly* result = malloc(sizeof(struct OrePoly));
  if (result == NULL)
  {
    printf("Could not allocate memory for the result of the multiplication.\n");
    exit(1);
  }
  result->degD1      = inp1->degD1 + inp2->degD1;
  result->degD2      = inp1->degD2 + inp2->degD2;
  result->ptrD1manip = inp1->ptrD1manip;
  result->ptrD2manip = inp1->ptrD2manip;
  result->coeffs     =
    malloc((result->degD1+1)*(result->degD2+1)*sizeof(struct GFModulus));
  if (result->coeffs == NULL)
  {
    printf("Could not allocate the coefficient vector for \
multiplication result\n");
    exit(1);
  }
  for (i = 0; i<(result->degD1+1)*(result->degD2+1); ++i)
    result->coeffs[i] = getZeroElemGF();
  if (result->degD1 + result->degD2 >= MIN_TOTAL_DEG_FOR_PARALLEL)
  {//condition fulfilled so that we go parallel
    int tempIter;
    omp_lock_t *resultLocks;
    resultLocks = (omp_lock_t*)malloc((inp1->degD2 + inp2->degD2 + 1)*sizeof(omp_lock_t));
    for (i=0; i<=inp1->degD2 + inp2->degD2;++i)
      omp_init_lock(resultLocks+i);
#pragma omp parallel for private(tempCoeff,m,i,j,k,l)
    for (tempIter = 0; tempIter<(inp1->degD1+1)*(inp1->degD2+1);++tempIter)
    {//iteration over the degrees of inp1
      i = tempIter/(inp1->degD1 +1);
      j = tempIter%(inp1->degD1 +1);
      for (k = 0; k<=inp2->degD2; ++k)
      {//iteration over degD2 of inp2
        omp_set_lock(resultLocks + (i+k));
        for (l = 0; l<=inp2->degD1; ++l)
        {//iteration over degD1 of inp2
          tempCoeff = inp2->coeffs[k*(inp2->degD1+1) + l];
#ifndef ORDERHOM2
          for (m = 0; m<i; ++m)
            tempCoeff = inp1->ptrD2manip(tempCoeff);
#else
          for (m = 0; m<i%ORDERHOM2; ++m)
            tempCoeff = inp1->ptrD2manip(tempCoeff);
#endif
#ifndef ORDERHOM1
          for (m = 0; m<j; ++m)
            tempCoeff = inp1->ptrD1manip(tempCoeff);
#else
          for (m = 0; m<j%ORDERHOM1; ++m)
            tempCoeff = inp1->ptrD1manip(tempCoeff);
#endif
          result->coeffs[(i+k)*(result->degD1+1) + j+l] =
            addGF(result->coeffs[(i+k)*(result->degD1+1) + j+l],
                  multGF((inp1->coeffs[i*(inp1->degD1+1)+j]),tempCoeff));
        }//iteration over degD1 of inp2
        omp_unset_lock(resultLocks + (i+k));
      }//iteration over degD2 of inp2
    }//iteration over the degrees of inp1
    for (i=0; i<=inp1->degD2 + inp2->degD2;++i)
      omp_destroy_lock(resultLocks+i);
    free(resultLocks);
  }//condition fulfilled so that we go parallel
  else
  {//non-parallel multiplication
    for (i = 0 ; i<=inp1->degD2; ++i)
    {//iteration over degD2 of inp1
      for (j = 0 ; j<=inp1->degD1; ++j)
      {//iteration over degD1 of inp1
        for (k = 0; k<=inp2->degD2; ++k)
        {//iteration over degD2 of inp2
          for (l = 0; l<=inp2->degD1; ++l)
          {//iteration over degD1 of inp2
            tempCoeff = inp2->coeffs[k*(inp2->degD1+1) + l];
#ifndef ORDERHOM2
            for (m = 0; m<i; ++m)
              tempCoeff = inp1->ptrD2manip(tempCoeff);
#else
            for (m = 0; m<i%ORDERHOM2; ++m)
              tempCoeff = inp1->ptrD2manip(tempCoeff);
#endif
#ifndef ORDERHOM1
            for (m = 0; m<j; ++m)
              tempCoeff = inp1->ptrD1manip(tempCoeff);
#else
            for (m = 0; m<j%ORDERHOM1; ++m)
              tempCoeff = inp1->ptrD1manip(tempCoeff);
#endif
            result->coeffs[(i+k)*(result->degD1+1) + j+l] =
              addGF(result->coeffs[(i+k)*(result->degD1+1) + j+l],
                    multGF((inp1->coeffs[i*(inp1->degD1+1)+j]),tempCoeff));
          }//iteration over degD1 of inp2
        }//iteration over degD2 of inp2
      }//iteration over degD1 of inp1
    }//iteration over degD2 of inp1
  }//non-parallel multiplication
  return result;
}//mult

struct OrePoly * getRandomOrePoly(int degD1, int degD2,
				  struct GFModulus
				  (*ptrD1manip)(struct GFModulus),
				  struct GFModulus
				  (*ptrD2manip)(struct GFModulus))
{//getRandomOrePoly
  int i;
  struct GFModulus * coeffs = malloc((degD1+1)*(degD2+1) *sizeof(struct GFModulus));
  for (i = 0; i<(degD1+1)*(degD2+1); ++i)
    coeffs[i] = getRandomGFElem();
  struct OrePoly* result = malloc(sizeof(struct OrePoly));
  if (result==NULL)
  {
    printf("Could not generate random polynomial, as no memory could \
be reserved.\n");
    exit(1);
  }
  result->degD1 = degD1;
  result->degD2 = degD2;
  result->coeffs = coeffs;
  result->ptrD1manip = ptrD1manip;
  result->ptrD2manip = ptrD2manip;
  return result;
}//getRandomOrePoly

struct OrePoly * generateRandomSecretKey(int deg, struct OrePoly* inp)
{//generateRandomSecretKey
  int i;

  struct OrePoly* result = malloc(sizeof(struct OrePoly));
  if (result == NULL)
  {
    printf("Cannot allocate memory for the result in \
generateRandomSecretKey.\n");
    exit(1);
  }
  result->degD1 = 0;
  result->degD2 = 0;
  struct GFModulus * coeffs = malloc(sizeof(struct GFModulus));
  if (coeffs == NULL)
  {
    printf("Could not reserve memory for one integer in \
scalarMult.\n");
    exit(1);
  }
  coeffs[0].coeffs[0]=(rand() % (MODULUS-1))+1;
  coeffs[0].coeffs[1]=0;
  coeffs[0].coeffs[2]=0;
  
  result->coeffs = coeffs;
  result->ptrD1manip = inp->ptrD1manip;
  result->ptrD2manip = inp->ptrD2manip;
  struct OrePoly * tempPPower=scalarMult(1,inp);//Instead of making a
						//lot of lines for copying
  struct OrePoly * temp1;
  struct OrePoly * temp2;
  int tempCoeff;
  for (i= 1 ; i<=deg; ++i)
  {
    if (i == deg)
      tempCoeff = (rand() % (MODULUS-1))+1;
    else
      tempCoeff = rand() % MODULUS;
    temp1 = scalarMult(tempCoeff,tempPPower); //c_i * P^i
    temp2 = add(result, temp1);//result + c_i*P^i
    free(temp1->coeffs);
    free(temp1);
    free(result->coeffs);
    free(result);
    result = temp2;
    temp1 = mult(tempPPower, inp);
    free(tempPPower->coeffs);
    free(tempPPower);
    tempPPower = temp1;
  }
  free(tempPPower->coeffs);
  free(tempPPower);
  return result;
}//generateRandomSecretKey

/* struct OrePoly** getAllPossibleOrePolysOfDegreeN(int n) */
/* {//getAllPossibleOrePolysOfDegreeN */
/*   if (n <0) */
/*     return NULL; */
/*   int expectedSize = 1; */
/*   int i; */
/*   for (i = 0; i<(n+1); ++i) */
/*     expectedSize*=NUMBEROFELEMENTSINGF; */
/*   struct GFModulus* allPossibleGFs = getAllPossibleGFElements(); */
/*   struct OrePoly** result = malloc(expectedSize*sizeOf(struct OrePoly*)); */
/*   if (n == 0) */
/*   {//termination case for the recursion */
/*     struct GFModulus *coeffs; */
/*     for (i = 0; i<NUMBEROFELEMENTSINGF; ++i) */
/*     { */
/*       result[i] = malloc(sizeof(struct OrePoly)); */
/*       result[i]->degD1 = 0; */
/*       result[i]->degD2 = 0; */
/*       result[i]->ptrD1manip = &Hom1; */
/*       result[i]->ptrD2manip = &Hom2; */
/*       coeffs = malloc(sizeof(GFModulus)); */
/*       coeffs[0] = allPossibleGFs[i]; */
/*       result[i]->coeffs = coeffs;       */
/*     } */
/*   }//termination case for the recursion */
/*   else */
/*   {//Recursive call case */
/*     struct GFModulus** resultTemp = getAllPossibleGFElements(n-1); */
/*     //TODO */
/*   return result; */
/* }//getAllPossibleOrePolysOfDegreeN */
