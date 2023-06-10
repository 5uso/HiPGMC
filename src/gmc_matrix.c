#include "gmc_matrix.h"

matrix new_matrix(int w, int h) {
    matrix m;
    m.w = w; m.h = h;
    m.data = malloc(w * h * sizeof(double));

    return m;
}

void free_matrix(matrix m) {
    free(m.data);
}

void print(matrix m) {
    for(int j = 0; j < m.h; j++) {
        for(int i = 0; i < m.w; i++) printf("%10.4lf", m.data[j * m.w + i]);
        printf("\n");
    }
}

sparse_matrix new_sparse(uint w, uint h) {
    sparse_matrix m;
    m.w = w; m.h = h;
    m.data = malloc(w * h * sizeof(sprs_val));

    return m;
}

void free_sparse(sparse_matrix m) {
    free(m.data);
}
