#include "funs.h"

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

matrix updateU(matrix q) { //Height of q is m
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

void updateF(matrix F, matrix U, double * ev, uint c) {
    for(int x = 0; x < U.w; x++) {
        double sum = 0.0d;
        for(int y = 0; y <= x; y++) {
            double t = -(U.data[y * U.w + x] + U.data[x * U.w + y]) / 2.0d;
            F.data[y * F.w + x] = t;
            sum += t;
        }
        F.data[x * F.w + x] -= sum;
    }

    //Eigenvalues go in ascending order inside w, eigenvectors are returned inside l. w is an array sized as a row of l.
    LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'U', F.w, F.data, F.w, ev); //TODO: Row major only outputs the upper triangular of eigenvectors?
}

int connectedComp(matrix m, int * y) {
    //TODO
    
}