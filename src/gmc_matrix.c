#include "gmc_matrix.h"

matrix new_matrix(int w, int h) {
    matrix m = {
        .w = w, .h = h,
        .data = malloc((long long) w * (long long) h * sizeof(double)),
    };
    return m;
}

void free_matrix(matrix m) {
    free(m.data);
}

void print(matrix m) {
    for(long long j = 0; j < m.h; j++) {
        for(long long i = 0; i < m.w; i++) printf("%10.4lf", m.data[j * m.w + i]);
        printf("\n");
    }
}

sparse_matrix new_sparse(int w, int h) {
    sparse_matrix m = {
        .w = w, .h = h,
        .data = malloc((long long) w * (long long) h * sizeof(sprs_val)),
    };
    return m;
}

void free_sparse(sparse_matrix m) {
    free(m.data);
}
