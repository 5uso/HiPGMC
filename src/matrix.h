#ifndef GMCMAT
#define GMCMAT

#include <stdlib.h>
#include <stdio.h>

typedef struct matrix {
    uint w; uint h;     //Width & Height
    double *  data;     //Matrix contents
} matrix;

matrix newMatrix(uint w, uint h);
void freeMatrix(matrix m);
double get(matrix m, uint i, uint j);
void set(matrix m, uint i, uint j, double val);
double * getRow(matrix m, uint j);
void print(matrix m);

#endif