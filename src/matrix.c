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

void print(matrix m) {
    for(int j = 0; j < m.h; j++) {
        for(int i = 0; i < m.w; i++) printf("%10.4lf", m.data[j * m.w + i]);
        printf("\n");
    }
}