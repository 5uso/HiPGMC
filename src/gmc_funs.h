#ifndef GMCFUNS
#define GMCFUNS

#include "gmc_scalapack.h"
#include "gmc_matrix.h"
#include "gmc_heap.h"
#include "gmc_sum.h"

#include <stdbool.h>
#include <lapacke.h>
#include <cblas.h>
#include <math.h>
#include <mpi.h>

#define EPS 2.2204460492503131e-16
#define BLOCK_SIZE 64

matrix sqr_dist(matrix m);
matrix update_u(matrix q);
matrix update_f(matrix F, matrix U, double * ev, int c, int rank, int blacs_row, int blacs_col, int blacs_height, int blacs_width, int blacs_ctx, MPI_Comm comm);
int connected_comp(bool * adj, int * y, int num);

#endif
