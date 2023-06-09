#ifndef GMCMAT
#define GMCMAT

#include <stdlib.h>
#include <stdio.h>

typedef struct matrix {
    uint w; uint h;     // Width & Height
    double *  data;     // Matrix contents
} matrix;

typedef struct sprs_val {
    uint i; double value;
} sprs_val;

typedef struct sparse_matrix {
    uint w; uint h;     // Width & Height (width is the number of non zero elements per row)
    sprs_val *data;     // Matrix contents
} sparse_matrix;

matrix new_matrix(uint w, uint h);
void free_matrix(matrix m);
void print(matrix m);
sparse_matrix new_sparse(uint w, uint h);
void free_sparse(sparse_matrix m);

#endif
