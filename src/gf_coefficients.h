#ifndef GF_COEFFICIENTS_H
#define GF_COEFFICIENTS_H

#define MODULUS (5)
#define DEGREEEXTENSION (3)
#define NUMBEROFELEMENTSINGF (125)

struct GFModulus
{//GFModulus
  int coeffs[DEGREEEXTENSION];
};//GFModulus

char* GFModulusToString(struct GFModulus);

void GFModulusToStdOut(struct GFModulus);

struct GFModulus addGF(struct GFModulus,struct GFModulus);

struct GFModulus minusGF(struct GFModulus, struct GFModulus);

struct GFModulus multGF(struct GFModulus, struct GFModulus);

int isZero_GF(struct GFModulus);

int isEqual_GF(struct GFModulus, struct GFModulus);

struct GFModulus getZeroElemGF(void);

struct GFModulus getIdentityElemGF(void);

struct GFModulus getMinusOneElemGF(void);

struct GFModulus getRandomGFElem(void);

struct GFModulus identityMap(struct GFModulus);

struct GFModulus Hom1(struct GFModulus);
#define ORDERHOM1 (3)

struct GFModulus Hom2(struct GFModulus);
#define ORDERHOM2 (3)

struct GFModulus scalarMultGF(int, struct GFModulus);

struct GFModulus* getAllPossibleGFElements(void);

#endif
