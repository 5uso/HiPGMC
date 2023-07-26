#include "gmc_scalapack.h"

static inline int _max(int a, int b) {
    return b > a ? b : a;
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

MPI_Datatype gmc_contiguous_long(MPI_Datatype type, long long count) {
    MPI_Datatype result;

    long long blocks = count / INT_MAX + 1;
    long long remainder = count % INT_MAX;

    MPI_Count element_size;
    MPI_Type_size_x(type, &element_size);
    
    int lengths[blocks]; MPI_Aint displs[blocks], d = 0; MPI_Datatype types[blocks];
    for(int i = 0; i < blocks; i++) {
        lengths[i] = INT_MAX;
        displs[i] = d;
        types[i] = type;
        d += INT_MAX * element_size;
    }
    lengths[blocks - 1] = remainder;

    // Construct a datatype with repetitions of the original type
    MPI_Type_create_struct(blocks, lengths, displs, types, &result);

    return result;
}

void gmc_distribute(int w, int h, double * a, double * b, int blacs_row, int blacs_col, int blacs_width, int blacs_height, int nb, int rank, MPI_Comm comm) {
    // Configure 2D block-cyclic distribution
    int numprocs = blacs_width * blacs_height, izero = 0;
    int global_dims[] = { h, w };
    int distr_modes[] = { MPI_DISTRIBUTE_CYCLIC, MPI_DISTRIBUTE_CYCLIC };
    int distr_block[] = { nb, nb };
    int procgr_dims[] = { blacs_width, blacs_height };

    // Dimensions of local matrix
    int mp = numroc_(&w, &nb, &blacs_row, &izero, &blacs_height);
    int nq = numroc_(&h, &nb, &blacs_col, &izero, &blacs_width);

    if(!rank) {
        #pragma omp parallel for
        for(int r = numprocs - 1; r >= 0; r--) {
            MPI_Datatype distr_arr; MPI_Request req;
            MPI_Type_create_darray(numprocs, r, 2, global_dims, distr_modes, distr_block, procgr_dims, MPI_ORDER_C, MPI_DOUBLE, &distr_arr);
            MPI_Type_commit(&distr_arr);

            // Send corresponding data to process
            if(!r) MPI_Isend(a, 1, distr_arr, r, 2711, comm, &req);
            else MPI_Send(a, 1, distr_arr, r, 2711, comm);

            MPI_Type_free(&distr_arr);
        }
    }

    // Receive from root process as contiguous array
    MPI_Status stat;
    MPI_Datatype cont = gmc_contiguous_long(MPI_DOUBLE, mp * (long long) nq);
    MPI_Type_commit(&cont);
    MPI_Recv(b, 1, cont, 0, 2711, comm, &stat);
    MPI_Type_free(&cont);
}

void gmc_collect(int w, int h, double * a, double * b, int blacs_row, int blacs_col, int blacs_width, int blacs_height, int nb, int rank, MPI_Comm comm) {
    // Configure 2D block-cyclic distribution
    int numprocs = blacs_width * blacs_height, izero = 0;
    int global_dims[] = { h, w };
    int distr_modes[] = { MPI_DISTRIBUTE_CYCLIC, MPI_DISTRIBUTE_CYCLIC };
    int distr_block[] = { nb, nb };
    int procgr_dims[] = { blacs_width, blacs_height };

    // Dimensions of local matrix
    int mp = numroc_(&w, &nb, &blacs_row, &izero, &blacs_height);
    int nq = numroc_(&h, &nb, &blacs_col, &izero, &blacs_width);

    // Create type to support potentially huge buffers
    MPI_Datatype cont = gmc_contiguous_long(MPI_DOUBLE, mp * (long long) nq);
    MPI_Type_commit(&cont);

    if(!rank) {
        MPI_Request req;
        MPI_Isend(a, 1, cont, 0, 2712, comm, &req);

        #pragma omp parallel for
        for(int r = 0; r < numprocs; r++) {
            MPI_Datatype distr_arr;
            MPI_Type_create_darray(numprocs, r, 2, global_dims, distr_modes, distr_block, procgr_dims, MPI_ORDER_C, MPI_DOUBLE, &distr_arr);
            MPI_Type_commit(&distr_arr);

            // Receive from each process and place into the global matrix
            MPI_Status stat;
            MPI_Recv(b, 1, distr_arr, r, 2712, comm, &stat);
        }

        return;
    }

    // Send data to root process as contiguous array
    MPI_Send(a, 1, cont, 0, 2712, comm);
    MPI_Type_free(&cont);
}
