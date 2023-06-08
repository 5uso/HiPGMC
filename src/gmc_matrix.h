#ifndef GMCMAT
#define GMCMAT

#include <stdlib.h>
#include <stdio.h>

typedef struct matrix {
    int w; int h;     // Width & Height
    double *  data;     // Matrix contents
} matrix;

matrix new_matrix(int w, int h);
void free_matrix(matrix m);
double get(matrix m, int i, int j);
void set(matrix m, int i, int j, double val);
double * get_row(matrix m, int j);
void print(matrix m);

#endif
