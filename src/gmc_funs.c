#include "gmc_funs.h"

inline int max(int a, int b) {
    return b > a ? b : a;
}

matrix sqr_dist(matrix m) {
    // Compute sum of squared columns vector
    double * ssc = malloc(m.w * sizeof(double));
    for(int i = 0; i < m.w; i++)
        ssc[i] = block_sum_col_sqr(m.data + i, m.h, m.w);

    // Compute multiplication by transpose (upper triangular only)
    matrix mt = new_matrix(m.w, m.w);
    cblas_dsyrk(CblasRowMajor, CblasUpper, CblasTrans, m.w, m.h, 1.0, m.data, m.w, 0.0, mt.data, m.w);

    // Compute final matrix
    matrix d = new_matrix(m.w, m.w);
    for(int i = 0; i < m.w; i++) {
        for(int j = 0; j < m.w; j++) {
            if(i == j) {
                d.data[j * d.w + i] = d.data[i * d.w + j] = 0.0;
                continue;
            }

            double mul = i < j ? mt.data[i * mt.w + j] : mt.data[j * mt.w + i];
            d.data[j * d.w + i] = d.data[i * d.w + j] = ssc[i] + ssc[j] - 2.0 * mul;
        }
    }

    free(ssc);
    free_matrix(mt);
    return d;
}

matrix update_u(matrix q) { // Height of q is m
    for(int j = 1; j < q.h; j++) {
        for(int i = 0; i < q.w; i++) q.data[i] += q.data[j * q.w + i];
    }

    double mean = 0.0;
    for(int i = 0; i < q.w; i++) mean += q.data[i];
    mean /= q.w;

    double vmin = INFINITY;
    for(int i = 0; i < q.w; i++) {
        q.data[i] = (q.data[i] - mean) / q.h + 1.0 / q.w;
        if(vmin > q.data[i]) vmin = q.data[i];
    }

    if(vmin >= 0.0) return q;

    double f = 1.0;
    double lambda_m = 0.0;
    double lambda_d = 0.0;
    for(int i = 0; i < q.w; i++) q.data[i] = -q.data[i];

    for(int ft = 1; ft <= 100 && (f < 0.0 ? -f : f) > 10e-10; ft++) {
        double sum = 0.0;
        int npos = 0;
        for(int i = 0; i < q.w; i++) {
            if((q.data[i] += lambda_d) > 0.0) {
                sum += q.data[i];
                npos++;
            } 
        }

        double g = (double)npos / q.w - 1.0 + EPS;
        f = sum / q.w - lambda_m;
        lambda_m += (lambda_d = -f / g);
    }

    for(int i = 0; i < q.w; i++) {
        if((q.data[i] = -q.data[i]) < 0.0) q.data[i] = 0.0;
    }

    q.h = 1;
    return q;
}

matrix update_f(matrix F, matrix U, double * ev, int c, int rank, int blacs_row, int blacs_col, int blacs_height, int blacs_width, int blacs_ctx,
                MPI_Comm comm, elpa_t handle, int * counts, int * displs) {
    int izero = 0, ione = 1, nb = BLOCK_SIZE, info;
    double done = 1.0, dzero = 0.0;
    arr_desc fd, flocald, eigvecd;

    int n = F.w;
    MPI_Bcast(&n, 1, MPI_INT, 0, comm);

    // Gather U into rank 0
    MPI_Gatherv(U.data, counts[rank], MPI_DOUBLE, F.data, counts, displs, MPI_DOUBLE, 0, comm);

    if(!rank) {
        for(int y = 0; y < F.w; y++) {
            for(int x = y; x < F.w; x++)
                F.data[y * F.w + x] = F.data[x * F.w + y] = (F.data[y * F.w + x] + F.data[x * F.w + y]) / -2.0;
                
            F.data[y * F.w + y] -= block_sum(F.data + y * F.w + y + 1, F.w - y - 1) + block_sum_col(F.data + y, y, F.w);
        }
    }
    
    // Find dimensions of distributed matrix and create local instances
    int mp = numroc_(&n, &nb, &blacs_row, &izero, &blacs_height);
    int nq = numroc_(&n, &nb, &blacs_col, &izero, &blacs_width);
    matrix f_local = new_matrix(mp, nq);
    int lld = max(1, numroc_(&n, &n, &blacs_row, &izero, &blacs_height));
    int lld_local = max(1, mp);

    // Distribute matrix by adding 0 to it
    // Original F has block size equal to total matrix size, since it's only on process 0
    // Distributed F is properly configured to be shared between processes
    descinit_(&fd, &n, &n, &n, &n, &izero, &izero, &blacs_ctx, &lld, &info);
    descinit_(&flocald, &n, &n, &nb, &nb, &izero, &izero, &blacs_ctx, &lld_local, &info);
    pdgeadd_("N", &n, &n, &done, F.data, &ione, &ione, &fd, &dzero, f_local.data, &ione, &ione, &flocald);

    // Call eigensolver. Eigenvalues are given in ascending order. We are responsible for freeing the returned buffers.
    // We are using row major upper triangular, but since fortran uses column major, we indicate lower triangular,
    // which is equivalent in symmetric matrices
    double *eigenvectors, *eigenvalues; int error;
    //gmc_pdsyevx('L', n, f_local.data, flocald, 1, c + 1, 0.0, &eigenvalues, &eigenvectors, &eigvecd);
    eigvecd = flocald;
    eigenvalues = malloc(n * sizeof(double)), eigenvectors = malloc(mp * nq * sizeof(double));
    elpa_eigenvectors(handle, f_local.data, eigenvalues, eigenvectors, &error);

    // Collect eigenvectors into process 0, set height to only work with c eigvecs, c+1 is only needed for eigval
    pdgeadd_("N", &n, &n, &done, eigenvectors, &ione, &ione, &eigvecd, &dzero, F.data, &ione, &ione, &fd);
    F.h = c;

    if(!rank) memcpy(ev, eigenvalues, (c + 1) * sizeof(double));

    free_matrix(f_local);
    free(eigenvectors);
    free(eigenvalues);
    return F;
}

int __find_comp(int * y, int i) {
    if(y[i] != i) y[i] = __find_comp(y, y[i]);
    return y[i];
}

void __merge_comp(int * y, int * ranks, int a, int b) {
    if(ranks[a] < ranks[b]) {
        y[a] = b;
        return;
    } 

    if(ranks[b] < ranks[a]) {
        y[b] = a;
        return;
    }

    y[b] = a;
    ranks[a]++;
}

int connected_comp(bool * adj, int * y, int num) {
    int * parents = malloc(num * sizeof(int));
    int * ranks = malloc(num * sizeof(int));
    memset(y, 0xFF, num * sizeof(int));

    for(int i = 0; i < num; i++) {
        parents[i] = i;
        ranks[i] = 0;
    }
    
    for(int j = 0; j < num; j++)
        for(int x = 0; x < j; x++)
            if(adj[j * num + x]) {
                int parentJ = __find_comp(parents, j);
                int parentX = __find_comp(parents, x);
                __merge_comp(parents, ranks, parentJ, parentX);
            } 

    int c = 0;
    for(int i = 0; i < num; i++) {
        int root = __find_comp(parents, i);
        if(y[root] == -1) y[root] = c++;
        y[i] = y[root];
    }

    free(parents);
    free(ranks);
    return c;
}
