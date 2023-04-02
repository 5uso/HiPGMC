#include "funs.h"

matrix sqr_dist(matrix m) {
    // Compute sum of squared columns vector
    double * ssc = malloc(m.w * sizeof(double));
    for(uint i = 0; i < m.w; i++) {
        double a = m.data[i];
        ssc[i] = a * a;
        for(uint j = 1; j < m.h; j++) {
            a = m.data[j * m.w + i];
            ssc[i] += a * a;
        }
    }

    // Compute multiplication by transpose (upper triangular only)
    matrix mt = new_matrix(m.w, m.w);
    cblas_dsyrk(CblasRowMajor, CblasUpper, CblasTrans, m.w, m.h, 1.0, m.data, m.w, 0.0, mt.data, m.w);

    // Compute final matrix
    matrix d = new_matrix(m.w, m.w);
    for(uint i = 0; i < m.w; i++) {
        for(uint j = 0; j < m.w; j++) {
            if(i == j) {
                d.data[j * d.w + i] = d.data[i * d.w + j] = 0.0d;
                continue;
            }

            double mul = i < j ? mt.data[i * mt.w + j] : mt.data[j * mt.w + i];
            d.data[j * d.w + i] = d.data[i * d.w + j] = ssc[i] + ssc[j] - 2.0d * mul;
        }
    }

    free(ssc);
    free_matrix(mt);
    return d;
}

matrix init_sig(matrix x, uint k) {
    matrix d = sqr_dist(x);

    for(int j = 0; j < x.w; j++) {
        heap h = new_heap(d.data + j * d.w, k + 2);
        for(int i = k + 2; i < x.w; i++) {
            if(d.data[j * d.w + i] < heap_max(h)) {
                *h.data[0] = 0.0;
                replace(&h, d.data + i + j * x.w);
            }
            else d.data[j * d.w + i] = 0.0d;
        }

        double a = heap_max(h);
        double b = a * k + EPS;
        for(int i = 1; i < k + 2; i++) {
            if(h.data[i] != h.min) b -= *h.data[i];
            else *h.data[i] = 0.0d;
        }

        for(int i = 0; i < k + 2; i++) {
            if(h.data[i] != h.min) *h.data[i] = (a - *h.data[i]) / b;
        }

        free_heap(h);
    }

    return d;
}

matrix update_u(matrix q) { // Height of q is m
    for(int j = 1; j < q.h; j++) {
        for(int i = 0; i < q.w; i++) q.data[i] += q.data[j * q.w + i];
    }

    double mean = 0.0d;
    for(int i = 0; i < q.w; i++) mean += q.data[i];
    mean /= q.w;

    double vmin = INFINITY;
    for(int i = 0; i < q.w; i++) {
        q.data[i] = (q.data[i] - mean) / q.h + 1.0d / q.w;
        if(vmin > q.data[i]) vmin = q.data[i];
    }

    if(vmin >= 0.0d) return q;

    double f = 1.0d;
    double lambda_m = 0.0d;
    double lambda_d = 0.0d;
    for(int i = 0; i < q.w; i++) q.data[i] = -q.data[i];

    for(int ft = 1; ft <= 100 && (f < 0.0d ? -f : f) > 10e-10; ft++) {
        double sum = 0.0d;
        int npos = 0;
        for(int i = 0; i < q.w; i++) {
            if((q.data[i] += lambda_d) > 0.0d) {
                sum += q.data[i];
                npos++;
            } 
        }

        double g = (double)npos / q.w - 1.0d + EPS;
        f = sum / q.w - lambda_m;
        lambda_m += (lambda_d = -f / g);
    }

    for(int i = 0; i < q.w; i++) {
        if((q.data[i] = -q.data[i]) < 0.0d) q.data[i] = 0.0d;
    }

    q.h = 1;
    return q;
}

matrix update_f(matrix F, matrix U, double * ev, uint c) {
    for(int y = 0; y < U.w; y++) {
        F.data[y * F.w + y] = -U.data[y * U.w + y];
        for(int i = 0; i < y; i++) F.data[y * F.w + y] -= F.data[y * U.w + i];
        for(int x = y + 1; x < U.w; x++) {
            double t = -(U.data[y * U.w + x] + U.data[x * U.w + y]) / 2.0d;
            F.data[x * F.w + y] = t;
            F.data[y * F.w + y] -= t;
        }
    }

    // Eigenvalues go in ascending order inside ev, eigenvectors are returned inside F
    /*LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'U', F.w, F.data, F.w, ev);
    F.h = c;
    return F;//*/

    int found_eigenvalue_n;
    int * ifail = malloc(sizeof(int) * c);
    double * eigenvectors = malloc(sizeof(double) * F.w * c);
    LAPACKE_dsyevx(LAPACK_COL_MAJOR, 'V', 'I', 'U', F.w, F.data, F.w, -1.0d, -1.0d, 1, c, 0.0d, &found_eigenvalue_n, ev, eigenvectors, F.w, ifail);
    F.h = c;
    memcpy(F.data, eigenvectors, sizeof(double) * F.w * F.h);
    free(ifail);
    free(eigenvectors);
    return F;
}

int __merge_comp(int * p, int c) {
    while(p[c] > c) c = p[c];
    return c;
}

int connected_comp(matrix m, int * y) {
    for(int i = 0; i < m.w; i++) y[i] = i;
    
    for(int j = 0; j < m.w; j++) {
        for(int x = j + 1; x < m.w; x++) {
            if(m.data[j * m.w + x] != 0.0d) y[__merge_comp(y, j)] = __merge_comp(y, x);
        }
    }

    int c = 0;
    for(int i = m.w - 1; i >= 0; i--) {
        if(i == y[i]) y[i] = c++;
        else y[i] = y[__merge_comp(y, y[i])];
    }

    return c;
}
