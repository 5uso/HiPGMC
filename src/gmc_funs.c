#include "gmc_funs.h"

matrix sqr_dist(matrix m, int rank, int blacs_row, int blacs_col, int blacs_height, int blacs_width, int blacs_ctx, MPI_Comm comm) {
    // Compute sum of squared columns vector
    double * ssc;
    if(!rank) {
        ssc = malloc(m.w * sizeof(double));
        #pragma omp parallel for
        for(int i = 0; i < m.w; i++)
            ssc[i] = block_sum_col_sqr(m.data + i, m.h, m.w);
    }

    // Sequential section, faster on some setups
    #ifdef SEQ_SQR
        if(rank) return m;
        matrix mt = new_matrix(m.w, m.w);
        cblas_dsyrk(CblasRowMajor, CblasUpper, CblasTrans, m.w, m.h, 1.0, m.data, m.w, 0.0, mt.data, m.w);
    #else
        // Distribute m
        MPI_Bcast(&m, sizeof(matrix), MPI_BYTE, 0, comm);
        int nb = BLOCK_SIZE, izero = 0, ione = 1, info;
        int mp = numroc_(&m.w, &nb, &blacs_row, &izero, &blacs_height);
        int nq = numroc_(&m.h, &nb, &blacs_col, &izero, &blacs_width);
        matrix m_local = new_matrix(mp, nq);
        gmc_distribute(m.w, m.h, m.data, m_local.data, blacs_row, blacs_col, blacs_width, blacs_height, nb, rank, comm);
        
        arr_desc mlocald, mtlocald;
        int lld_local = mp > 1 ? mp : 1;
        descinit_(&mlocald, &m.w, &m.h, &nb, &nb, &izero, &izero, &blacs_ctx, &lld_local, &info);
        nq = numroc_(&m.w, &nb, &blacs_col, &izero, &blacs_width);
        descinit_(&mtlocald, &m.w, &m.w, &nb, &nb, &izero, &izero, &blacs_ctx, &lld_local, &info);
        
        // Compute multiplication by transpose (upper triangular only)
        matrix mt_local = new_matrix(mp, nq);
        double done = 1.0, dzero = 0.0;
        pdsyrk_("L", "N", &m.w, &m.h, &done, m_local.data, &ione, &ione, &mlocald, &dzero, mt_local.data, &ione, &ione, &mtlocald);
        free_matrix(m_local);

        // Collect mt
        matrix mt;
        if(!rank) mt = new_matrix(m.w, m.w);
        gmc_collect(m.w, m.w, mt_local.data, mt.data, blacs_row, blacs_col, blacs_width, blacs_height, nb, rank, comm);
        free_matrix(mt_local);

        // Workers can return here
        if(rank) return mt;
    #endif

    // Compute final matrix
    matrix d = new_matrix(m.w, m.w);
    #pragma omp parallel for
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
    int izero = 0, nb = BLOCK_SIZE, info;
    arr_desc fd, flocald;

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

    // Distribute matrix by adding 0 to it
    // Original F has block size equal to total matrix size, since it's only on process 0
    // Distributed F is properly configured to be shared between processes
    int lld_local = mp > 1 ? mp : 1;
    descinit_(&fd, &n, &n, &n, &n, &izero, &izero, &blacs_ctx, &n, &info);
    descinit_(&flocald, &n, &n, &nb, &nb, &izero, &izero, &blacs_ctx, &lld_local, &info);
    gmc_distribute(n, n, F.data, f_local.data, blacs_row, blacs_col, blacs_width, blacs_height, nb, rank, comm);

    // Call eigensolver. Eigenvalues are given in ascending order. We are responsible for freeing the returned buffers.
    // We are using row major upper triangular, but since fortran uses column major, we indicate lower triangular,
    // which is equivalent in symmetric matrices
    double *eigenvectors, *eigenvalues;
    #ifdef ELPA_API_VER
        int error;
        eigenvalues = malloc(n * sizeof(double)), eigenvectors = malloc(mp * nq * sizeof(double));
        elpa_eigenvectors(handle, f_local.data, eigenvalues, eigenvectors, &error);
    #else
        arr_desc eigvecd;
        gmc_pdsyevx('L', n, f_local.data, flocald, 1, c + 1, 0.0, &eigenvalues, &eigenvectors, &eigvecd);
    #endif

    // Collect eigenvectors into process 0, set height to only work with c eigvecs, c+1 is only needed for eigval
    gmc_collect(n, n, eigenvectors, F.data, blacs_row, blacs_col, blacs_width, blacs_height, nb, rank, comm);
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
