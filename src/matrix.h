#ifndef GMCMAT
#define GMCMAT

#include <stdlib.h>
#include <stdio.h>

typedef struct matrix {
    uint w; uint h;     //Width & Height
    double *  data;     //Matrix contents
} matrix;

matrix new_matrix(uint w, uint h);
void free_matrix(matrix m);
double get(matrix m, uint i, uint j);
void set(matrix m, uint i, uint j, double val);
double * get_row(matrix m, uint j);
void print(matrix m);

#endif