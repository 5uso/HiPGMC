#ifndef GMC
#define GMC

#include "funs.h"

typedef struct gmc_result {
    matrix y;
    matrix U;
    matrix S0;
    matrix F;
    matrix evs;
} gmc_result;

gmc_result gmc(matrix * X, uint c, uint lambda, bool normalize);

#endif