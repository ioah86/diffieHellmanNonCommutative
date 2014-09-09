/*
  Running the Diffie-Hellman protocol as stated in the paper 
  "A Diffie-Hellman-like Key Exchange Protocol Based on Multivariate Ore Polynomials"
  (c) Reinhold Burger and Albert Heinle
 */

#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include "ore_algebra.c"

#define LDEG 25,25
#define QDEG 3,2
#define PDEG 2,3
#define PPDEG 10
#define QQDEG 10

int main()
{//main
//  int i;
  srand(time(NULL));

  /* struct OrePoly p; */
  /* struct OrePoly q; */
  /* p.degD1 = 2; */
  /* p.degD2 = 3; */
  /* q.degD1 = 3; */
  /* q.degD2 = 2; */
  /* p.ptrD1manip = &Hom1; */
  /* p.ptrD2manip = &Hom2; */
  /* q.ptrD1manip = &Hom1; */
  /* q.ptrD2manip = &Hom2; */
  
  /* int p00[3] = {2,1,1}; */
  /* int p10[3] = {4,3,2}; */
  /* int p20[3] = {1,3,3}; */
  /* int p01[3] = {3,0,4}; */
  /* int p11[3] = {3,2,1}; */
  /* int p21[3] = {3,2,4}; */
  /* int p02[3] = {0,4,2}; */
  /* int p12[3] = {3,3,3}; */
  /* int p22[3] = {3,2,0}; */
  /* int p03[3] = {4,3,3}; */
  /* int p13[3] = {0,0,4}; */
  /* int p23[3] = {3,1,4}; */

  /* int q00[3] = {1,4,2}; */
  /* int q10[3] = {4,4,0}; */
  /* int q20[3] = {0,3,2}; */
  /* int q30[3] = {1,3,1}; */
  /* int q01[3] = {2,3,1}; */
  /* int q11[3] = {0,4,1}; */
  /* int q21[3] = {0,2,0}; */
  /* int q31[3] = {0,3,3}; */
  /* int q02[3] = {0,3,0}; */
  /* int q12[3] = {4,2,3}; */
  /* int q22[3] = {1,0,3}; */
  /* int q32[3] = {3,4,2}; */


  /* struct GFModulus coeffsp[12]; */
  /* coeffsp[0].coeffs[0] = p00[0]; */
  /* coeffsp[0].coeffs[1] = p00[1]; */
  /* coeffsp[0].coeffs[2] = p00[2]; */
  /* coeffsp[1].coeffs[0] = p10[0]; */
  /* coeffsp[1].coeffs[1] = p10[1]; */
  /* coeffsp[1].coeffs[2] = p10[2]; */
  /* coeffsp[2].coeffs[0] = p20[0]; */
  /* coeffsp[2].coeffs[1] = p20[1]; */
  /* coeffsp[2].coeffs[2] = p20[2]; */

  /* coeffsp[3].coeffs[0] = p01[0]; */
  /* coeffsp[3].coeffs[1] = p01[1]; */
  /* coeffsp[3].coeffs[2] = p01[2]; */
  /* coeffsp[4].coeffs[0] = p11[0]; */
  /* coeffsp[4].coeffs[1] = p11[1]; */
  /* coeffsp[4].coeffs[2] = p11[2]; */
  /* coeffsp[5].coeffs[0] = p21[0]; */
  /* coeffsp[5].coeffs[1] = p21[1]; */
  /* coeffsp[5].coeffs[2] = p21[2]; */

  /* coeffsp[6].coeffs[0] = p02[0]; */
  /* coeffsp[6].coeffs[1] = p02[1]; */
  /* coeffsp[6].coeffs[2] = p02[2]; */
  /* coeffsp[7].coeffs[0] = p12[0]; */
  /* coeffsp[7].coeffs[1] = p12[1]; */
  /* coeffsp[7].coeffs[2] = p12[2]; */
  /* coeffsp[8].coeffs[0] = p22[0]; */
  /* coeffsp[8].coeffs[1] = p22[1]; */
  /* coeffsp[8].coeffs[2] = p22[2]; */

  /* coeffsp[9].coeffs[0] = p03[0]; */
  /* coeffsp[9].coeffs[1] = p03[1]; */
  /* coeffsp[9].coeffs[2] = p03[2]; */
  /* coeffsp[10].coeffs[0] = p13[0]; */
  /* coeffsp[10].coeffs[1] = p13[1]; */
  /* coeffsp[10].coeffs[2] = p13[2]; */
  /* coeffsp[11].coeffs[0] = p23[0]; */
  /* coeffsp[11].coeffs[1] = p23[1]; */
  /* coeffsp[11].coeffs[2] = p23[2]; */

  /* struct GFModulus coeffsq[12]; */

  /* coeffsq[0].coeffs[0] = q00[0]; */
  /* coeffsq[0].coeffs[1] = q00[1]; */
  /* coeffsq[0].coeffs[2] = q00[2]; */
  /* coeffsq[1].coeffs[0] = q10[0]; */
  /* coeffsq[1].coeffs[1] = q10[1]; */
  /* coeffsq[1].coeffs[2] = q10[2]; */
  /* coeffsq[2].coeffs[0] = q20[0]; */
  /* coeffsq[2].coeffs[1] = q20[1]; */
  /* coeffsq[2].coeffs[2] = q20[2]; */
  /* coeffsq[3].coeffs[0] = q30[0]; */
  /* coeffsq[3].coeffs[1] = q30[1]; */
  /* coeffsq[3].coeffs[2] = q30[2]; */

  /* coeffsq[4].coeffs[0] = q01[0]; */
  /* coeffsq[4].coeffs[1] = q01[1]; */
  /* coeffsq[4].coeffs[2] = q01[2]; */
  /* coeffsq[5].coeffs[0] = q11[0]; */
  /* coeffsq[5].coeffs[1] = q11[1]; */
  /* coeffsq[5].coeffs[2] = q11[2]; */
  /* coeffsq[6].coeffs[0] = q21[0]; */
  /* coeffsq[6].coeffs[1] = q21[1]; */
  /* coeffsq[6].coeffs[2] = q21[2]; */
  /* coeffsq[7].coeffs[0] = q31[0]; */
  /* coeffsq[7].coeffs[1] = q31[1]; */
  /* coeffsq[7].coeffs[2] = q31[2]; */

  /* coeffsq[8].coeffs[0] = q02[0]; */
  /* coeffsq[8].coeffs[1] = q02[1]; */
  /* coeffsq[8].coeffs[2] = q02[2]; */
  /* coeffsq[9].coeffs[0] = q12[0]; */
  /* coeffsq[9].coeffs[1] = q12[1]; */
  /* coeffsq[9].coeffs[2] = q12[2]; */
  /* coeffsq[10].coeffs[0] = q22[0]; */
  /* coeffsq[10].coeffs[1] = q22[1]; */
  /* coeffsq[10].coeffs[2] = q22[2]; */
  /* coeffsq[11].coeffs[0] = q32[0]; */
  /* coeffsq[11].coeffs[1] = q32[1]; */
  /* coeffsq[11].coeffs[2] = q32[2]; */

  /* p.coeffs = coeffsp; */
  /* q.coeffs = coeffsq; */

  /* char* pInStr = OrePolyToString(&p); */
  /* char* qInStr = OrePolyToString(&q); */
  /* printf("P: %s\n", pInStr); */
  /* printf("Q: %s\n", qInStr); */
  /* free(pInStr); */
  /* free(qInStr); */
  
  /* int i; */

  /* struct OrePoly *tempPPow1 = mult(&p,&p); */
  /* struct OrePoly *tempPPow2; */
  /* struct OrePoly *tempQPow1 = mult(&q,&q); */
  /* struct OrePoly *tempQPow2; */
  
  /* for (i = 2; i<=50; i++) */
  /* { */
  /*   pInStr = OrePolyToString(tempPPow1); */
  /*   qInStr = OrePolyToString(tempQPow1); */
  /*   printf("P^%i: %s\n",i,pInStr); */
  /*   printf("Q^%i: %s\n",i,qInStr); */
  /*   free(pInStr); */
  /*   free(qInStr); */
  /*   tempPPow2 = mult(tempPPow1,&p); */
  /*   tempQPow2 = mult(tempQPow1,&q); */
  /*   free(tempPPow1); */
  /*   free(tempQPow1); */
  /*   tempPPow1 = tempPPow2; */
  /*   tempQPow1 = tempQPow2; */
  /* } */

  //public Part
  struct OrePoly* L =
    getRandomOrePoly(LDEG,&Hom1,&Hom2);
  struct OrePoly* P =
    getRandomOrePoly(PDEG,&Hom1,&Hom2);
  struct OrePoly* Q =
    getRandomOrePoly(QDEG,&Hom1,&Hom2);
  char* LString = OrePolyToString(L);
  printf("L: %s\n", LString);
  free(LString);
  char* PString = OrePolyToString(P);
  printf("P: %s\n", PString);
  free(PString);
  char* QString = OrePolyToString(Q);
  printf("Q: %s\n", QString);
  free(QString);

  //Private Part: Alice
  printf("Calculating P_A...\n");
  struct OrePoly* P_A = generateRandomSecretKey(PPDEG,P);
  char *P_AToString = OrePolyToString(P_A);
  printf("Done. Result: P_A = %s\n", P_AToString);
  free(P_AToString);
  printf("Calculating Q_A...\n");
  struct OrePoly* Q_A = generateRandomSecretKey(QQDEG,Q);
  char *Q_AToString = OrePolyToString(Q_A);
  printf("Done. Result: Q_A = %s\n", Q_AToString);
  free(Q_AToString);

  //Private Part: Bob
  printf("Calculating P_B...\n");
  struct OrePoly* P_B = generateRandomSecretKey(PPDEG,P);
  char *P_BToString = OrePolyToString(P_B);
  printf("Done. Result: P_B = %s\n", P_BToString);
  free(P_BToString);
  printf("Calculating Q_B...\n");
  struct OrePoly* Q_B = generateRandomSecretKey(QQDEG,Q);
  char *Q_BToString = OrePolyToString(Q_B);
  printf("Done. Result: Q_B = %s\n", Q_BToString);
  free(Q_BToString);

  //Sending Stuff
  //Alice to Bob
  printf("Calculating Message Alice to Bob:\n");
  struct OrePoly* temp = mult(P_A, L);
  struct OrePoly* AToB = mult(temp, Q_A);
  free(temp->coeffs);
  free(temp);
  char * AToBToString = OrePolyToString(AToB);
  printf("Done. Result is: %s\n", AToBToString);
  free(AToBToString);
  
  //Bob to Alice
  printf("Calculating Message Bob to Alice:\n");
  temp = mult(P_B,L);
  struct OrePoly* BToA = mult(temp, Q_B);
  free(temp->coeffs);
  free(temp);
  char * BToAToString = OrePolyToString(BToA);
  printf("Done. Result is: %s\n", BToAToString);
  free(BToAToString);
  
  printf("Calculating the secret key:\n");
  //Secret Key computation
  temp = mult(P_B, AToB);
  struct OrePoly* sKey = mult(temp, Q_B);
  free(temp->coeffs);
  free(temp);
  char * sKeyToString = OrePolyToString(sKey);
  printf("Done. The secret key is: %s\n", sKeyToString);
  free(sKeyToString);

  free(L->coeffs);
  free(L);
  free(P->coeffs);
  free(P);
  free(Q->coeffs);
  free(Q);
  free(P_A->coeffs);
  free(P_A);
  free(P_B->coeffs);
  free(P_B);
  free(Q_A->coeffs);
  free(Q_A);
  free(Q_B->coeffs);
  free(Q_B);
  free(AToB->coeffs);
  free(AToB);
  free(BToA->coeffs);
  free(BToA);
  free(sKey->coeffs);
  free(sKey);
  return 0;
}//main

