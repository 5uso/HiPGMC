#ifndef GMCFUNS
#define GMCFUNS

#include "gmc_matrix.h"
#include "gmc_heap.h"
#include "gmc_sum.h"

#include <stdbool.h>
#include <lapacke.h>
#include <cblas.h>
#include <math.h>

#define EPS 2.2204460492503131e-16

matrix sqr_dist(matrix m);
matrix init_sig(matrix x, uint k);
matrix update_u(matrix q);
matrix update_f(matrix F, matrix U, double * ev, uint c);
int connected_comp(bool * adj, int * y, int num);

#endif
