#ifndef GMCSCALAPACK
#define GMCSCALAPACK

#include <stdlib.h>
#include <math.h>

typedef struct arr_desc {
    int dtype, ctxt, m, n, mb, nb, rsrc, csrc, lld;
} arr_desc;

extern void Cblacs_pinfo(int*, int*);
extern void Cblacs_get(int, int, int*);
extern void Cblacs_gridinit(int*, const char*, int, int);
extern void Cblacs_gridinfo(int, int*, int*, int*, int*);
extern void Cblacs_pcoord(int, int, int*, int*);
extern void Cblacs_gridexit(int);
extern void Cblacs_barrier(int, const char*);
extern void Cdgerv2d(int, int, int, double*, int, int, int);
extern void Cdgesd2d(int, int, int, double*, int, int, int);

extern int numroc_(int*, int*, int*, int*, int*);
extern void descinit_(arr_desc*, int*, int*, int*, int*, int*, int*, int*, int*, int*);
extern void pdgeadd_(char*, int*, int*, double*, double*, int*, int*, arr_desc*, double*, double*, int*, int*, arr_desc*);
extern void pdsyevx_(char * jobz, char * range, char * uplo, int * n, double * a, int * ia, int * ja, arr_desc * desca, double * vl, double * vu, int * il, int * iu, double * abstol, int * m, int * nz, double * w, double * orfac, double * z, int * iz, int * jz, arr_desc * descz, double * work, int * lwork, int * iwork, int * liwork, int * ifail, int * iclustr, double * gap, int * info);
extern double pdlange_(char * norm, int * m, int * n, double * a, int * ia, int * ja, arr_desc * desca, double * work);

int gmc_pdsyevx(char uplo, int n, double * a, arr_desc desca, int il, int iu, double abstol, double ** w, double ** z, arr_desc * descz);

#endif
