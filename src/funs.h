#ifndef GMCFUNS
#define GMCFUNS

#include "matrix.h"
#include "heap.h"
#include <math.h>
#include <cblas.h>
#include <lapacke.h>

#define EPS 2.2204e-16

matrix sqr_dist(matrix m);
matrix init_sig(matrix x, uint k);
matrix update_u(matrix q);
matrix update_f(matrix F, matrix U, double * ev, uint c);
int connected_comp(matrix m, int * y);

#endif
