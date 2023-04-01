#ifndef GMC
#define GMC

#include "funs.h"
#include <math.h>
#include <string.h>
#include <stdbool.h>

typedef struct gmc_result {
    matrix U;
    matrix * S0;
    matrix F;
    matrix evs;
    int * y;
    int n;
    int m;
    int cluster_num;
    int iterations;
    double lambda;
} gmc_result;

gmc_result gmc(matrix * X, uint m, uint c, double lambda, bool normalize);
void free_gmc_result(gmc_result r);

#define NITER 20
#define ZR 10e-11
#define PN 15
#define IS_LOCAL 1
#define INLINE_GMC_INTERNALS 1

#ifdef INLINE_GMC_INTERNALS
    #define GMC_INTERNAL static inline
#else
    #define GMC_INTERNAL static
#endif

#endif
