#ifndef GMCMAT
#define GMCMAT

#include <stdlib.h>
#include <stdio.h>

typedef struct matrix {
    int w; int h; // Width & Height
    double *data; // Matrix contents
} matrix;

typedef struct sprs_val {
    int i; double value;
} sprs_val;

typedef struct sparse_matrix {
    int w; int h;   // Width & Height (width is the number of non zero elements per row)
    sprs_val *data; // Matrix contents
} sparse_matrix;

matrix new_matrix(int w, int h);
void free_matrix(matrix m);
void print(matrix m);
sparse_matrix new_sparse(int w, int h);
void free_sparse(sparse_matrix m);

#endif
