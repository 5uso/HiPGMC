#include "gmc_scale.h"

static inline int _max(int a, int b) {
    return b > a ? b : a;
}

static inline long long _minl(long long a, long long b) {
    return b < a ? b : a;
}

int gmc_pdsyevx(char uplo, int n, double * a, arr_desc desca, int il, int iu, double abstol, double ** w, double ** z, arr_desc * descz) {
    // Blacs info
    int ctx = desca.ctxt, nb = desca.nb, blacs_height, blacs_width, blacs_row, blacs_col;
    Cblacs_gridinfo(ctx, &blacs_height, &blacs_width, &blacs_row, &blacs_col);

    // Call variables
    int izero = 0, ione = 1;
    double dzero = 0.0;

    int m, nz, info;

    // Prepare Z
    int zm = n, zn = n;
    int z_mp = numroc_(&zm, &nb, &blacs_row, &izero, &blacs_height);
    int z_nq = numroc_(&zn, &nb, &blacs_col, &izero, &blacs_width);
    *z = malloc(z_mp * z_nq * sizeof(double));
    int z_lld_distr = _max(1, z_mp);
    descinit_(descz, &zm, &zn, &nb, &nb, &izero, &izero, &ctx, &z_lld_distr, &info);

    // Other buffers
    int lwork = -1, liwork = -1;
    int * ifail = malloc(n * sizeof(int));
    int * iwork = malloc(100 * sizeof(int));
    int * icluster = malloc(2 * blacs_height * blacs_width * sizeof(int));
    double * work = malloc(100 * sizeof(double));
    double * gap = malloc(blacs_height * blacs_width * sizeof(double));
    *w = malloc(zn * sizeof(double));

    // Call to query workspace
    pdsyevx_("V", "I", &uplo, &n, a, &ione, &ione, &desca, &dzero, &dzero, &il, &iu, &abstol, &m, &nz, *w, &dzero, *z, &ione, &ione, descz, work, &lwork, iwork, &liwork, ifail, icluster, gap, &info);
    if(info) return info;
    
    // Prepare workspace
    int nn = _max(_max(n, nb), 2);
    int np = numroc_(&nn, &nb, &izero, &izero, &blacs_height);
    int mq = numroc_(&nn, &nb, &izero, &izero, &blacs_width);
    int nnp = _max(_max(n, 4), blacs_height * blacs_width + 1);
    lwork = _max((int) fabs(work[0]), 5 * n + _max(5 * nn, np * mq + 2 * nb * nb) + ((n - 1) / (blacs_height * blacs_width) + 1) * nn);
    liwork = _max(iwork[0], 6 * nnp);
    work = realloc(work, lwork * sizeof(double));
    iwork = realloc(iwork, liwork * sizeof(int));

    // Finally, do the thing
    pdsyevx_("V", "I", &uplo, &n, a, &ione, &ione, &desca, &dzero, &dzero, &il, &iu, &abstol, &m, &nz, *w, &dzero, *z, &ione, &ione, descz, work, &lwork, iwork, &liwork, ifail, icluster, gap, &info);

    // Cleanup
    free(ifail);
    free(iwork);
    free(icluster);
    free(work);
    free(gap);

    return info;
}

static inline void _copy_cyclic(int w, int h, int nb, int blacs_width, int blacs_height, int blacs_col, int blacs_row, double * from,
                                double * to, long long max_bytes, long long * y, int * y_blk, long long * x, int * x_blk) {
    long long count = 0, wb = w * sizeof(double), nbb = nb * sizeof(double);
    char * a =  (char *) from, * b = (char *) to;
    for(; *y < h; *y += nb * (blacs_width - 1)) {
        for(; *y_blk < nb && *y < h; (*y_blk)++, (*y)++) {
            for(; *x < wb; *x += nbb * (blacs_height - 1)) {
                while(*x_blk < nbb && *x < wb) {
                    b[count++] = a[*y * wb + (*x)++]; (*x_blk)++;
                    if(count >= max_bytes) {
                        // If we fill the buffer, leave the loop with the proper values to do next iteration
                        if(*x_blk >= nbb || *x >= wb) { *x_blk = 0; *x += nbb * (blacs_height - 1); }
                        if(*x >= wb) { *x = (long long) (nb * blacs_row) * sizeof(double); (*y_blk)++, (*y)++; }
                        if(*y_blk >= nb || *y >= h) { *y_blk = 0;  *y += nb * (blacs_width - 1); }
                        return;
                    }
                }
                *x_blk = 0;
            }
            *x = (long long) (nb * blacs_row) * sizeof(double);
        }
        *y_blk = 0;
    }
}

static inline void _fill_cyclic(int w, int h, int nb, int blacs_width, int blacs_height, int blacs_col, int blacs_row, double * from,
                                double * to, long long max_bytes, long long * y, int * y_blk, long long * x, int * x_blk) {
    long long count = 0, wb = w * sizeof(double), nbb = nb * sizeof(double);
    char * a =  (char *) from, * b = (char *) to;
    for(; *y < h; *y += nb * (blacs_width - 1)) {
        for(; *y_blk < nb && *y < h; (*y_blk)++, (*y)++) {
            for(; *x < wb; *x += nbb * (blacs_height - 1)) {
                while(*x_blk < nbb && *x < wb) {
                    b[*y * wb + (*x)++] = a[count++]; (*x_blk)++;
                    if(count >= max_bytes) {
                        // If we fill the buffer, leave the loop with the proper values to do next iteration
                        if(*x_blk >= nbb || *x >= wb) { *x_blk = 0; *x += nbb * (blacs_height - 1); }
                        if(*x >= wb) { *x = (long long) (nb * blacs_row) * sizeof(double); (*y_blk)++, (*y)++; }
                        if(*y_blk >= nb || *y >= h) { *y_blk = 0;  *y += nb * (blacs_width - 1); }
                        return;
                    }
                }
                *x_blk = 0;
            }
            *x = (long long) (nb * blacs_row) * sizeof(double);
        }
        *y_blk = 0;
    }
}

void gmc_distribute(int w, int h, double * a, double * b, int rank, int blacs_width, int blacs_height, int nb, MPI_Comm comm) {
    int izero = 0, numprocs = blacs_width * blacs_height;
    MPI_Bcast(&w, 1, MPI_INT, 0, comm);
    MPI_Bcast(&h, 1, MPI_INT, 0, comm);

    if(!rank) {
        #pragma omp parallel for
        for(int r = 0; r < numprocs; r++) {
            // Dimensions of r's local matrix
            int blacs_col = r / blacs_height;
            int blacs_row = r % blacs_height;
            long long mp = numroc_(&w, &nb, &blacs_row, &izero, &blacs_height);
            long long nq = numroc_(&h, &nb, &blacs_col, &izero, &blacs_width);
            long long numbytes = mp * nq * sizeof(double);

            // Set up block-cyclic distribution start
            long long y = nb * blacs_col, x = (long long) (nb * blacs_row) * sizeof(double);
            int y_blk = 0, x_blk = 0;
            if(!r) {
                // Self: copy into local buffer
                _copy_cyclic(w, h, nb, blacs_width, blacs_height, blacs_col, blacs_row, a, b, LLONG_MAX, &y, &y_blk, &x, &x_blk);
                continue;
            }

            // Send to process: split the message in chunks to control max size
            double * buf = malloc(_minl(numbytes, MAX_MPI_MSG_BYTES));
            for(long long pos = 0; pos < numbytes; pos += MAX_MPI_MSG_BYTES) {
                _copy_cyclic(w, h, nb, blacs_width, blacs_height, blacs_col, blacs_row, a, buf, MAX_MPI_MSG_BYTES, &y, &y_blk, &x, &x_blk);
                long long amt = _minl(numbytes - pos, MAX_MPI_MSG_BYTES);
                MPI_Send(buf, (int) amt, MPI_BYTE, r, 2711, comm);
            }
            free(buf);
        }

        return;
    }

    // Dimensions of local matrix
    int blacs_col = rank / blacs_height;
    int blacs_row = rank % blacs_height;
    long long mp = numroc_(&w, &nb, &blacs_row, &izero, &blacs_height);
    long long nq = numroc_(&h, &nb, &blacs_col, &izero, &blacs_width);
    long long numbytes = mp * nq * sizeof(double);

    // Receive matching chunks
    char * recvbuf = (char *) b;
    for(long long pos = 0; pos < numbytes; pos += MAX_MPI_MSG_BYTES) {
        long long amt = _minl(numbytes - pos, MAX_MPI_MSG_BYTES);
        MPI_Recv(recvbuf + pos, (int) amt, MPI_BYTE, 0, 2711, comm, MPI_STATUS_IGNORE);
    }
}

void gmc_collect(int w, int h, double * a, double * b, int rank, int blacs_width, int blacs_height, int nb, MPI_Comm comm) {
    int izero = 0, numprocs = blacs_width * blacs_height;
    MPI_Bcast(&w, 1, MPI_INT, 0, comm);
    MPI_Bcast(&h, 1, MPI_INT, 0, comm);

    if(!rank) {
        #pragma omp parallel for
        for(int r = 0; r < numprocs; r++) {
            // Dimensions of r's local matrix
            int blacs_col = r / blacs_height;
            int blacs_row = r % blacs_height;
            long long mp = numroc_(&w, &nb, &blacs_row, &izero, &blacs_height);
            long long nq = numroc_(&h, &nb, &blacs_col, &izero, &blacs_width);
            long long numbytes = mp * nq * sizeof(double);

            // Set up block-cyclic distribution start
            long long y = nb * blacs_col, x = (long long) (nb * blacs_row) * sizeof(double);
            int y_blk = 0, x_blk = 0;
            if(!r) {
                // Self: place values from local buffer into global matrix
                _fill_cyclic(w, h, nb, blacs_width, blacs_height, blacs_col, blacs_row, a, b, LLONG_MAX, &y, &y_blk, &x, &x_blk);
                continue;
            }

            // Receive from process: receive matching chunks and place them into global matrix
            double * buf = malloc(_minl(numbytes, MAX_MPI_MSG_BYTES));
            for(long long pos = 0; pos < numbytes; pos += MAX_MPI_MSG_BYTES) {
                long long amt = _minl(numbytes - pos, MAX_MPI_MSG_BYTES);
                MPI_Recv(buf, (int) amt, MPI_BYTE, r, 2712, comm, MPI_STATUS_IGNORE);
                _fill_cyclic(w, h, nb, blacs_width, blacs_height, blacs_col, blacs_row, buf, b, MAX_MPI_MSG_BYTES, &y, &y_blk, &x, &x_blk);
            }
            free(buf);
        }

        return;
    }

    // Dimensions of local matrix
    int blacs_col = rank / blacs_height;
    int blacs_row = rank % blacs_height;
    long long mp = numroc_(&w, &nb, &blacs_row, &izero, &blacs_height);
    long long nq = numroc_(&h, &nb, &blacs_col, &izero, &blacs_width);
    long long numbytes = mp * nq * sizeof(double);

    // Send message in chunks to ensure we don't exceed max message size
    char * sendbuf = (char *) a;
    for(long long pos = 0; pos < numbytes; pos += MAX_MPI_MSG_BYTES) {
        long long amt = _minl(numbytes - pos, MAX_MPI_MSG_BYTES);
        MPI_Send(sendbuf + pos, (int) amt, MPI_BYTE, 0, 2712, comm);
    }
}

void gmc_scatter(long long w, int h, void * a, void * b, int rank, int numprocs, MPI_Comm comm) {
    MPI_Bcast(&w, 1, MPI_LONG, 0, comm);
    MPI_Bcast(&h, 1, MPI_INT, 0, comm);

    if(!rank) {
        #pragma omp parallel for
        for(int r = 0; r < numprocs; r++) {
            // Rows assigned to process
            long long numrows = h / numprocs + (r < h % numprocs);
            long long numbytes = numrows * w;
            if(!r) {
                // Self: copy directly into local
                memcpy(b, a, numbytes);
                continue;
            }

            // Send to process: split message in chunks
            void * offset = a + ((long long) (h / numprocs) * (long long) r + _minl(h % numprocs, r)) * w;
            for(long long pos = 0; pos < numbytes; pos += MAX_MPI_MSG_BYTES) {
                long long amt = _minl(numbytes - pos, MAX_MPI_MSG_BYTES);
                MPI_Send(offset + pos, (int) amt, MPI_BYTE, r, 2713, comm);
            }
        }

        return;
    }

    // Rows assigned to process
    long long numrows = h / numprocs + (rank < h % numprocs);
    long long numbytes = numrows * w;

    // Receive matching chunks
    for(long long pos = 0; pos < numbytes; pos += MAX_MPI_MSG_BYTES) {
        long long amt = _minl(numbytes - pos, MAX_MPI_MSG_BYTES);
        MPI_Recv(b + pos, (int) amt, MPI_BYTE, 0, 2713, comm, MPI_STATUS_IGNORE);
    }
}

void gmc_gather(long long w, int h, void * a, void * b, int rank, int numprocs, MPI_Comm comm) {
    MPI_Bcast(&w, 1, MPI_LONG, 0, comm);
    MPI_Bcast(&h, 1, MPI_INT, 0, comm);

    if(!rank) {
        #pragma omp parallel for
        for(int r = 0; r < numprocs; r++) {
            // Rows assigned to process
            long long numrows = h / numprocs + (r < h % numprocs);
            long long numbytes = numrows * w;
            if(!r) {
                // Self: copy directly from local
                memcpy(b, a, numbytes);
                continue;
            }

            // Receive matching chunks from process
            void * offset = b + ((long long) (h / numprocs) * (long long) r + _minl(h % numprocs, r)) * w;
            for(long long pos = 0; pos < numbytes; pos += MAX_MPI_MSG_BYTES) {
                long long amt = _minl(numbytes - pos, MAX_MPI_MSG_BYTES);
                MPI_Recv(offset + pos, (int) amt, MPI_BYTE, r, 2714, comm, MPI_STATUS_IGNORE);
            }
        }

        return;
    }

    // Rows assigned to process
    long long numrows = h / numprocs + (rank < h % numprocs);
    long long numbytes = numrows * w;

    // Send message in chunks
    for(long long pos = 0; pos < numbytes; pos += MAX_MPI_MSG_BYTES) {
        long long amt = _minl(numbytes - pos, MAX_MPI_MSG_BYTES);
        MPI_Send(a + pos, (int) amt, MPI_BYTE, 0, 2714, comm);
    }
}
