#include "matrix.h"

matrix newMatrix(uint w, uint h) {
    matrix m;
    m.w = w; m.h = h;
    m.data = malloc(w * h * sizeof(double));

    return m;
}

void freeMatrix(matrix m) {
    free(m.data);
}

double get(matrix m, uint i, uint j) { // i = x axis, j = y axis
    return m.data[j * m.w + i];
}

void set(matrix m, uint i, uint j, double val) {
    m.data[j * m.w + i] = val;
}

double * getRow(matrix m, uint j) {
    return m.data + j * m.w;
}

void print(matrix m) {
    for(int j = 0; j < m.h; j++) {
        for(int i = 0; i < m.w; i++) printf("%lf ", get(m, i, j));
        printf("\n");
    }
}