/*
  Implementation of the Diffie-Hellman protocol -- Non-commutative version

  (c) Reinhold Burger and Albert Heinle
 */

#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>

#define MODULUS (5)
#define DEGREEEXTENSION (3)

//////////////////////////////////////////////////
// CAUTION: The Following functions and definitions have to be made
// mathematically correct; i.e. a person who modifies this part of the
// code has to know what he/she is doing.
// We assume the minimal polynomial of a to be
// x^3 + 3*x + 3
//////////////////////////////////////////////////
struct GFModulus
{//GFModulus
  int coeffs[DEGREEEXTENSION];
};//GFModulus

char* GFModulusToString(struct GFModulus inp)
{
  int coeffInStringSize = MODULUS/10+1;
  char* result = malloc(DEGREEEXTENSION*coeffInStringSize + 13);
  sprintf(result,"(%i + %ia + %ia^2)",inp.coeffs[0],inp.coeffs[1],
		  inp.coeffs[2]);
  return result;
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

char* OrePolyToString(struct OrePoly* inp)
{//printOrePoly
  int i; int j;//iteration variables
  int coeffInStringSize = DEGREEEXTENSION*(MODULUS/10+1) + 13;
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
}//printOrePoly

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

struct OrePoly *getIdentityElemOrePoly()
{//getIdentityElemOrePoly
  struct OrePoly *result = malloc(sizeof(struct OrePoly));
  result->degD1 = 0;
  result->degD2 = 0;
  result->ptrD1manip = &Hom1;
  result->ptrD2manip = &Hom2;
  struct GFModulus *coeffs = malloc(sizeof(struct GFModulus));
  coeffs[0] = getIdentityElemGF();
  result->coeffs = coeffs;
  return result;
}//getIdentityElemOrePoly

struct OrePoly *getZeroElemOrePoly()
{//getZeroElemOrePoly
  struct OrePoly *result = malloc(sizeof(struct OrePoly));
  result->degD1 = 0;
  result->degD2 = 0;
  result->ptrD1manip = &Hom1;
  result->ptrD2manip = &Hom2;
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
    struct GFModulus zeroElem;
    zeroElem.coeffs[0]= 0;
    zeroElem.coeffs[1]= 0;
    zeroElem.coeffs[2]= 0;//TODO make it dependend on DEGREEEXTENSION
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
  //memset(result->coeffs,0,(result->degD1+1)*(result->degD2+1)*sizeof(int));
  for (i = 0 ; i<=inp1->degD2; ++i)
  {
    for (j = 0 ; j<=inp1->degD1; ++j)
    {
      for (k = 0; k<=inp2->degD2; ++k)
      {
	for (l = 0; l<=inp2->degD1; ++l)
	{
//	  printf("i: %i , j: %i , k: %i , l:%i\n",i,j,k,l);
//	  printf("position: %i\n",(i+k)*(result->degD1+1) + j+l);
	  tempCoeff = inp2->coeffs[k*(inp2->degD1+1) + l];
	  for (m = 0; m<i; ++m)
	    tempCoeff = inp1->ptrD2manip(tempCoeff);
	  for (m = 0; m<j; ++m)
	    tempCoeff = inp1->ptrD1manip(tempCoeff);
//	  printf("We are adding %i * %i to that position.\n",
//		 (inp1->coeffs[i*(inp1->degD1+1)+j]),tempCoeff );
//	  printf("initial value at that position: %i\n",
//		 result->coeffs[(i+k)*(result->degD1+1) + j+l]);
	  result->coeffs[(i+k)*(result->degD1+1) + j+l] =
	    addGF(result->coeffs[(i+k)*(result->degD1+1) + j+l],
		  multGF((inp1->coeffs[i*(inp1->degD1+1)+j]),tempCoeff));
//	  printf("Value afterwards: %i\n",
//		 result->coeffs[(i+k)*(result->degD1+1) + j+l]);
	}
      }
    }
  }
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
  coeffs[0].coeffs[0]=1;
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
    tempCoeff = rand() % MODULUS;
    if (i == deg)
      if (tempCoeff == 0)//making sure degree is correct
	tempCoeff=1;
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

