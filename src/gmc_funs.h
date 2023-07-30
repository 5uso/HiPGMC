#ifndef GMCFUNS
#define GMCFUNS

#include "gmc_scalapack.h"
#include "gmc_matrix.h"
#include "gmc_heap.h"
#include "gmc_sum.h"

#include <elpa/elpa.h>
#include <stdbool.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>

#define EPS 2.2204460492503131e-16
#define BLOCK_SIZE 64
#define ELPA_API_VER 20221109
//#define ELPA_GPU

// Run square matrix multiplication on a single node, which may be faster on some systems
//#define SEQ_SQR

matrix sqr_dist(matrix m, int rank, int blacs_row, int blacs_col, int blacs_height, int blacs_width, int blacs_ctx, MPI_Comm comm);
matrix update_u(matrix q);
matrix update_f(matrix F, double * ev, int c, int rank, int blacs_row, int blacs_col, int blacs_height, int blacs_width, int blacs_ctx, MPI_Comm comm, elpa_t handle);
int connected_comp(bool * adj, int * y, int num);

#endif
