/*
  Running the Diffie-Hellman protocol as stated in the paper 
  "A Diffie-Hellman-like Key Exchange Protocol Based on Multivariate Ore Polynomials"
  (c) Reinhold Burger and Albert Heinle
 */

#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include "ore_algebra.h"

#define LDEG 5,5
#define QDEG 3,2
#define PDEG 2,3
#define PPDEG 2
#define QQDEG 2

int main()
{//main
//  int i;
  srand(time(NULL));

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

