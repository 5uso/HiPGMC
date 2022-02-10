#include "funs.h"

#include <stdio.h>

matrix sqrDist(matrix m) {
    //Compute sum of squared columns vector
    double * ssc = malloc(m.w * sizeof(double));
    for(uint i = 0; i < m.w; i++) {
        double a = get(m, i, 0);
        ssc[i] = a * a;
        for(uint j = 1; j < m.h; j++) {
            a = get(m, i, j);
            ssc[i] += a * a;
        }
    }

    //Compute multiplication by transpose (upper triangular only)
    matrix mt = newMatrix(m.w, m.w);
    cblas_dsyrk(CblasRowMajor, CblasUpper, CblasTrans, m.w, m.h, 1.0, m.data, m.w, 0.0, mt.data, m.w);

    //Compute final matrix
    matrix d = newMatrix(m.w, m.w);
    for(uint i = 0; i < m.w; i++) {
        for(uint j = 0; j < m.w; j++) {
            if(i == j) {
                set(d, i, j, 0.0d);
                continue;
            }

            double mul = i < j ? get(mt, j, i) : get(mt, i, j);
            set(d, i, j, ssc[i] + ssc[j] - 2.0d * mul);
        }
    }

    free(ssc);
    freeMatrix(mt);
    return d;
}

matrix initSIG(matrix x, uint k) {
    matrix d = sqrDist(x);

    for(int j = 0; j < x.w; j++) {
        heap h = newHeap(getRow(d, j), k + 2);
        for(int i = k + 2; i < x.w; i++) {
            if(get(d, i, j) < heapMax(h)) replace(&h, d.data + i + j * x.w);
            else set(d, i, j, 0.0d);
        }

        double a = heapMax(h);
        double b = a * k + EPS;
        for(int i = 1; i < k + 2; i++) {
            if(h.data[i] != h.min) b -= *h.data[i];
            else *h.data[i] = 0.0d;
        }

        for(int i = 0; i < k + 2; i++) {
            if(h.data[i] != h.min) *h.data[i] = (a - *h.data[i]) / b;
        }

        freeHeap(h);
    }

    return d;
}

matrix updateU(matrix q, uint m) {
    for(int j = 0; j < q.h; j++) {
        for(int i = 0; i < q.w; i++) {
            m.data[i] += m.data[j * m.w + i];
        }
    }

    double vmin = 0.0d;
    if(vmin >= 0.0d) return q;

    double f = 1.0d;
    double lambda_m = 0.0d;
    for(int ft = 1; ft <= 100 && (f < 0.0d ? -f : f) < 10e-10; ft++) {

    }
    return q;
}

matrix eig(matrix l, uint c) {
    matrix w = newMatrix(l.w, 1);
    LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'U', l.w, l.data, l.w, w.data); //TODO: Row major only outputs the upper triangular of eigenvectors?
    return w;
}