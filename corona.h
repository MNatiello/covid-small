/*
 * New corona.h for version 11. Minor changes in sizes and variable names
 * 2020-10-20
 */
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>

#define max(a,b)  ((a) > (b) ? (a) : (b))   /* Estos no se usan */
#define min(a,b)  ((a) < (b) ? (a) : (b))

#include "../dSFMT-src-2.2.3/dSFMT.h"
uint32_t seed=17;    /*Should be unsigned int. No point in using -1 or any negative*/
dsfmt_t dsfmt;
static const int dsfmt_mexp = DSFMT_MEXP;

#define Maxdays              730          /* 2 years */
#define MaxStat                5000       /* Max number of repetitions */
#define CONT                       9          /* contagion states */
#define SIZE         1+2*(CONT+1)   /* SIZE=S,T[CONT],U[CONT],TR,UR */
#define Events      1+6*CONT          /* extern CONT*{infectionT infectionU promotionT promotionU remotionT remotionU } */

/*Rate parameters*/
double          ext,betaT,betaU,detT,detU,remT[CONT],remU[CONT];
unsigned int delayT,delayU;
double          slow,eps;
double          vir[CONT],pr[CONT];
/*end rate parameters*/

unsigned int X[SIZE]; /* subpopulations */
unsigned int Cases;   /* acumulador de casos detectados */
double W[Events];     /* Event-rates */
int wipe,Tmin;     /* option to eliminate short-lived contagion cycles */
