#ifndef _GLOBAL_H
#define _GLOBAL_H

#include "chromosome.h"
#include "myrand.h"

extern bool SHOW_POPULATION;
extern bool SHOW_REPLACEMENT;
extern bool SHOW_SELECTION_INDEX;
extern double R1;
extern double R2;
extern double Bratio;
extern int NUM_CUTPOINTS;

extern bool isRTR;






extern void outputErrMsg (const char *errMsg);
extern MyRand myRand;
extern int numOfReplacement;
#endif
