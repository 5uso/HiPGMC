#ifndef GMC
#define GMC

#include "funs.h"
#include <string.h>

typedef struct gmc_result {
    matrix y;
    matrix U;
    matrix S0;
    matrix F;
    matrix evs;
} gmc_result;

gmc_result gmc(matrix * X, uint m, uint c, double lambda, bool normalize);

#define NITER 20
#define ZR 10e-11
#define PN 15
#define IS_LOCAL 1

#endif