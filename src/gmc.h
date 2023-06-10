#ifndef GMC
#define GMC

#include "gmc_scalapack.h"
#include "gmc_matrix.h"
#include "gmc_funs.h"
#include "gmc_heap.h"
#include "gmc_sum.h"

#include <elpa/elpa.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

typedef struct gmc_result {
    matrix U;
    sparse_matrix * S0;
    matrix F;
    matrix evs;
    int * y;
    int n;
    int m;
    int cluster_num;
    int iterations;
    double lambda;
} gmc_result;

gmc_result gmc(matrix * X, int m, int c, double lambda, bool normalize, MPI_Comm in_comm, int in_context);
void free_gmc_result(gmc_result r);

#define NITER 20
#define ZR 10e-11
#define PN 15
#define IS_LOCAL
#define INLINE_GMC_INTERNALS
//#define PRINT_GMC_STEPS

#ifdef INLINE_GMC_INTERNALS
    #define GMC_INTERNAL static inline
#else
    #define GMC_INTERNAL __attribute__((noinline))
#endif

#ifdef PRINT_GMC_STEPS
    #define GMC_STEP(x) if(!rank) (x)
#else
    #define GMC_STEP(x)
#endif

#endif
