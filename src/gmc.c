#include "gmc.h"

static int numprocs, rank;
static MPI_Comm comm;

static int blacs_row, blacs_col, blacs_height, blacs_width;
static int blacs_ctx;

static elpa_t handle;

static int * pattern_cnts;
static int * counts;
static int * displs;
static MPI_Datatype gmc_row_type;

#ifdef PRINT_GMC_STEPS
    #include <sys/time.h>
    struct timeval begin, curr;
#endif

GMC_INTERNAL void __gmc_normalize(matrix * X, int m, int num) {
    for(int v = 0; v < m; v++) {
        long long h = X[v].h, w = X[v].w;
        #pragma omp parallel for
        for(long long x = 0; x < w; x++) {
            double mean = 0.0;
            for(long long y = 0; y < h; y++) mean += X[v].data[y * w + x];
            mean /= h;

            double std = 0.0;
            for(long long y = 0; y < h; y++) {
                double dev = X[v].data[y * w + x] - mean;
                std += dev * dev;
            }
            std /= h - 1;
            std = sqrt(std);
            if(std == 0) std = EPS;

            for(long long y = 0; y < h; y++) X[v].data[y * w + x] = (X[v].data[y * w + x] - mean) / std;
        }
    }
}

GMC_INTERNAL void __gmc_init_s0(matrix * X, int m, int num, sparse_matrix * S0, matrix * ed, double * sums) {
    for(long long v = 0; v < m; v++) {
        matrix ted; // Rank 0 computes square distance matrix and scatters
        ted = sqr_dist(rank ? ted : X[v], rank, blacs_row, blacs_col, blacs_height, blacs_width, blacs_ctx, comm);
        matrix local_ted = new_matrix(num, pattern_cnts[rank]);
        MPI_Scatterv(ted.data, counts, displs, gmc_row_type, local_ted.data, counts[rank], gmc_row_type, 0, comm);
        if(!rank) free_matrix(ted);

        S0[v] = new_sparse(PN + 1, local_ted.h);
        ed[v] = new_matrix(PN + 1, local_ted.h);

        int s = displs[rank]; // Start pattern for this process
        #pragma omp parallel for
        for(long long y = 0; y < pattern_cnts[rank]; y++) {
            local_ted.data[y * num + s + y] = INFINITY;
            
            heap h = new_heap(local_ted.data + y * num, PN + 1);
            for(long long x = PN + 1; x < num; x++)
                if(local_ted.data[y * num + x] < heap_max(h))
                    replace(&h, local_ted.data + y * num + x);
                
            sums[v * pattern_cnts[rank] + y] = block_sum_ptr(h.data + 1, PN, 0);
            double denominator = *h.data[0] * PN - sums[v * pattern_cnts[rank] + y] + EPS;

            for(long long i = 0; i < PN + 1; i++) {
                sprs_val val = {
                    .i = h.data[i] - (local_ted.data + y * num),
                    .value = ((*h.data[0] - *h.data[i]) / denominator) * (i > 0),
                };
                S0[v].data[y * (PN + 1) + i] = val;
                ed[v].data[y * (PN + 1) + i] = *h.data[i];
            }

            free_heap(h);
        }

        free_matrix(local_ted);
    }
}

GMC_INTERNAL matrix __gmc_init_u(sparse_matrix * S0, int m, int num) {
    matrix U = new_matrix(num, pattern_cnts[rank]);

    // U starts as average of SIG matrices
    memset(U.data, 0x00, (long long) num * (long long) pattern_cnts[rank] * sizeof(double));

    #pragma omp parallel for
    for(long long y = 0; y < pattern_cnts[rank]; y++) {
        double sum = 0.0;
        for(long long i = 0; i < PN + 1; i++)
            for(int v = 0; v < m; v++) {
                sprs_val val = S0[v].data[y * (PN + 1) + i];
                double t = val.value / m;
                U.data[y * num + val.i] += t;
                sum += t;
            }
        for(long long x = 0; x < num; x++) U.data[y * num + x] /= sum;
    }

    return U;
}

GMC_INTERNAL void __gmc_update_s0(sparse_matrix * S0, matrix U, matrix w, int m, int num, matrix * ed, double * sums) {
    for(long long v = 0; v < m; v++) {
        double weight = w.data[v] * 2.0;

        #pragma omp parallel for
        for(long long y = 0; y < pattern_cnts[rank]; y++) {
            double max = ed[v].data[(PN + 1) * y];
            double maxU = U.data[y * num + S0[v].data[(PN + 1) * y].i];

            double sumU = 0.0;
            for(long long i = 1; i < PN + 1; i++) {
                long long x = S0[v].data[y * (PN + 1) + i].i;
                sumU += U.data[y * num + x];
            }

            double numerator = max - weight * maxU;
            double denominator = PN * max - sums[v * pattern_cnts[rank] + y] + weight * (sumU - PN * maxU) + EPS;

            for(long long i = 0; i < PN + 1; i++) {
                long long x = S0[v].data[y * (PN + 1) + i].i;
                double r = (numerator - ed[v].data[(PN + 1) * y + i] + weight * U.data[y * num + x]) / denominator;
                S0[v].data[y * (PN + 1) + i].value = r * (r > 0.0);
            }
        }
    }
}

GMC_INTERNAL void __gmc_update_w(sparse_matrix * S0, matrix U, matrix w, int m, int num) {
    matrix US = new_matrix(num, pattern_cnts[rank]);
    matrix US_global; if(!rank) US_global = new_matrix(num, num);

    // Set up distributed array descriptor
    arr_desc desca_local; int nb = BLOCK_SIZE, izero = 0, ione = 1, info;
    int mp = numroc_(&num, &nb, &blacs_row, &izero, &blacs_height);
    int nq = numroc_(&num, &nb, &blacs_col, &izero, &blacs_width);
    int lld_local = mp > 1 ? mp : 1;
    matrix US_local = new_matrix(mp, nq);
    descinit_(&desca_local, &num, &num, &nb, &nb, &izero, &izero, &blacs_ctx, &lld_local, &info);

    for(int v = 0; v < m; v++) {
        memcpy(US.data, U.data, (long long) num * (long long) pattern_cnts[rank] * sizeof(double));
        
        #pragma omp parallel for
        for(long long y = 0; y < pattern_cnts[rank]; y++)
            for(long long i = 0; i < PN + 1; i++) {
                sprs_val val = S0[v].data[y * (PN + 1) + i];
                long long x = val.i;
                US.data[y * num + x] -= val.value;
            }

        // Redistribute matrix from row to block cyclic
        MPI_Gatherv(US.data, counts[rank], gmc_row_type, US_global.data, counts, displs, gmc_row_type, 0, comm);
        gmc_distribute(num, num, US_global.data, US_local.data, blacs_row, blacs_col, blacs_width, blacs_height, nb, rank, comm);

        // Compute frobenius norm in parallel
        double distUS = pdlange_("F", &num, &num, US_local.data, &ione, &ione, &desca_local, NULL);
        w.data[v] = 0.5 / (distUS + EPS);
    }

    free_matrix(US); free_matrix(US_local);
    if(!rank) free_matrix(US_global);
}

GMC_INTERNAL void __gmc_update_u(sparse_matrix * S0, matrix U, matrix w, matrix * F, int m, int num, double * lambda) {
    matrix dist; // Rank 0 computes square distance matrix and scatters
    dist = sqr_dist(*F, rank, blacs_row, blacs_col, blacs_height, blacs_width, blacs_ctx, comm); // F is transposed, since LAPACK returns it in column major
    matrix local_dist = new_matrix(num, pattern_cnts[rank]);
    MPI_Scatterv(dist.data, counts, displs, gmc_row_type, local_dist.data, counts[rank], gmc_row_type, 0, comm);
    if(!rank) free_matrix(dist);

    #pragma omp parallel for
    for(long long y = 0; y < pattern_cnts[rank]; y++) {
        int qw = 0;
        int * idx = malloc((long long) num * sizeof(int));

        #ifdef IS_LOCAL
            memset(idx, 0x00, (long long) num * sizeof(int));
            for(int v = 0; v < m; v++) {
                for(long long i = 0; i < PN + 1; i++) {
                    sprs_val val = S0[v].data[y * (PN + 1) + i];
                    long long x = val.i;
                    qw -= idx[x];
                    idx[x] |= (val.value > 0);
                    qw += idx[x];
                }
            }
        #else
            memset(idx, 0x01, num * sizeof(int));
            qw = num;
        #endif

        matrix q = new_matrix(qw, m);
        for(long long x = 0, i = 0; x < num; x++) {
            if(idx[x]) {
                q.data[i] = *lambda * local_dist.data[y * num + x] / (double) m * -0.5;
                i++;
                idx[x] = i;
            } 
        }

        for(long long v = m - 1; v >= 0; v--)
            for(long long i = 0; i < qw; i++)
                q.data[v * qw + i] = q.data[i] / w.data[v];

        for(long long v = m - 1; v >= 0; v--) {
            for(long long i = 0; i < PN + 1; i++) {
                sprs_val val = S0[v].data[y * (PN + 1) + i];
                long long x = val.i;
                if(idx[x]) q.data[v * qw + idx[x] - 1] += val.value;
            }
        }

        q = update_u(q);
        for(long long x = 0, i = 0; x < num; x++) {
            if(idx[x]) {
                U.data[y * num + x] = q.data[i];
                i++;
            } else U.data[y * num + x] = 0.0;
        }

        free_matrix(q);
        free(idx);
    }

    free_matrix(local_dist);
}

GMC_INTERNAL bool __gmc_main_loop(int it, sparse_matrix * S0, matrix U, matrix w, matrix * F, matrix * F_old, matrix evs, int m, int c, int num,
                                  matrix * ed, double * sums, double * lambda) {
    bool end_iteration = false;
    double * ev = NULL;

    // Update S0
    GMC_STEP("update S0", it);
    __gmc_update_s0(S0, U, w, m, num, ed, sums);

    // Update w
    GMC_STEP("update w", it);
    __gmc_update_w(S0, U, w, m, num);

    // Update U
    GMC_STEP("update U", it);
    __gmc_update_u(S0, U, w, F, m, num, lambda);

    // Update matrix of eigenvectors (F), as well as eigenvalues
    GMC_STEP("update F", it);
    matrix temp = *F_old;
    *F_old = *F;
    *F = temp;
    ev = evs.data + (long long)(c + 1) * (long long) it;

    MPI_Gatherv(U.data, counts[rank], gmc_row_type, F->data, counts, displs, gmc_row_type, 0, comm);
    *F = update_f(*F, ev, c, rank, blacs_row, blacs_col, blacs_height, blacs_width, blacs_ctx, comm, handle);

    // Update lambda
    GMC_STEP("update lambda", it);
    if(!rank) {
        double fn = block_sum(ev, c);
        if(fn > ZR) {
            *lambda *= 2.0;
        } else if(fn + ev[c] < ZR) {
            *lambda /= 2.0;
            temp = *F_old;
            *F_old = *F;
            *F = temp;
        } else {
            evs.h = it + 2;
            end_iteration = true;
        }
    }
    
    struct { bool end; double l; } linfo = {
        .end = end_iteration,
        .l = *lambda,
    };

    MPI_Bcast(&linfo, sizeof(linfo), MPI_BYTE, 0, comm);
    end_iteration = linfo.end; *lambda = linfo.l;
    return end_iteration;
}

elpa_t configure_elpa(int n, int nev) {
    int error, izero = 0, nb = BLOCK_SIZE;
    if(elpa_init(ELPA_API_VER) != ELPA_OK) {
        fprintf(stderr, "{ %d, %d } Error: ELPA API version not supported\n", blacs_row, blacs_col);
        exit(1);
    }

    elpa_t handle = elpa_allocate(&error);
    elpa_set(handle, "na", n, &error);
    elpa_set(handle, "nev", nev, &error);
    elpa_set(handle, "nblk", nb, &error);
    elpa_set(handle, "mpi_comm_parent", MPI_Comm_c2f(comm), &error);
    elpa_set(handle, "process_row", blacs_row, &error);
    elpa_set(handle, "process_col", blacs_col, &error);
    elpa_set(handle, "blacs_context", blacs_ctx, &error);

    arr_desc descriptor;
    int mp = numroc_(&n, &nb, &blacs_row, &izero, &blacs_height);
    int nq = numroc_(&n, &nb, &blacs_col, &izero, &blacs_width);
    int lld = mp > 1 ? mp : 1;
    descinit_(&descriptor, &n, &n, &nb, &nb, &izero, &izero, &blacs_ctx, &lld, &error);
    if(error) {
        fprintf(stderr, "{ %d, %d } Error: Invalid distribution, descinit arg %d\n", blacs_row, blacs_col, -error);
        exit(1);
    }

    elpa_set(handle, "local_nrows", mp, &error);
    elpa_set(handle, "local_ncols", nq, &error);

    error = elpa_setup(handle);
    if(error != ELPA_OK) {
        fprintf(stderr, "{ %d, %d } Error: Elpa setup returned code %d\n", blacs_row, blacs_col, error);
        exit(1);
    }

    elpa_set(handle, "solver", ELPA_SOLVER_2STAGE, &error);
    if(!rank && error != ELPA_OK) fprintf(stderr, "Warning: Couldn't set ELPA 2stage solver\n");

    #ifdef ELPA_GPU
        elpa_set(handle, "nvidia-gpu", 1, &error);
        if(error == ELPA_OK) {
            elpa_set(handle, "real_kernel", ELPA_2STAGE_REAL_NVIDIA_GPU, &error);
            if(!rank && error != ELPA_OK) fprintf(stderr, "Warning: Couldn't set ELPA gpu kernel\n");
        } else if(!rank) {
            fprintf(stderr, "Warning: ELPA GPU acceleration not supported\n");
        }
    #else
        elpa_set(handle, "real_kernel", ELPA_2STAGE_REAL_AVX512_BLOCK2, &error);
        if(!rank && error != ELPA_OK) printf("Warning: Couldn't set ELPA AVX512 kernel\n");

        char * num_threads_env = getenv("OMP_NUM_THREADS");
        if(num_threads_env) {
            int thread_num = (int) strtol(num_threads_env, NULL, 10);
            if(thread_num) {
                elpa_set(handle, "omp_threads", thread_num, &error);
                if(!rank && error != ELPA_OK) fprintf(stderr, "Warning: Couldn't set ELPA OMP threads to %d\n", thread_num);
            }
        }
    #endif

    return handle;
}

gmc_result gmc(matrix * X, int m, int c, double lambda, bool normalize, MPI_Comm in_comm, int in_context) {
    #ifdef PRINT_GMC_STEPS
        gettimeofday(&begin, 0);
    #endif

    sparse_matrix *S0 = NULL; double *sums = NULL;   
    matrix *ed = NULL, U, F, F_old, evs, w;

    // Get MPI info
    comm = in_comm;
    MPI_Comm_size(comm, &numprocs);
    MPI_Comm_rank(comm, &rank);

    // Get BLACS info
    blacs_ctx = in_context;
    Cblacs_gridinfo(blacs_ctx, &blacs_height, &blacs_width, &blacs_row, &blacs_col);

    // Broadcast parameters to all processes in the group
    struct { int m, c, n; double l; bool norm; } params = {
        .m = m,
        .c = c,
        .n = rank ? 0 : X[0].w,
        .l = lambda,
        .norm = normalize,
    };

    MPI_Bcast(&params, sizeof(params), MPI_BYTE, 0, comm);
    m = params.m, c = params.c, lambda = params.l, normalize = params.norm;
    int num = params.n;

    #ifdef ELPA_API_VER
        handle = configure_elpa(num, c + 1);
    #endif

    if(!rank) {
        // Normalize data
        GMC_STEP("Init: normalize");
        if(normalize) __gmc_normalize(X, m, num);
    }

    // Determine how many patterns correspond to each process
    pattern_cnts = malloc((long long) numprocs * sizeof(int));
    counts = malloc((long long) numprocs * sizeof(int));
    displs = malloc((long long) numprocs * sizeof(int));
    int displacement = 0;
    for(int i = 0; i < numprocs; i++) {
        pattern_cnts[i] = num / numprocs + (i < num % numprocs);
        counts[i] = pattern_cnts[i];
        displs[i] = displacement;
        displacement += counts[i];
    }

    gmc_row_type = gmc_contiguous_long(MPI_DOUBLE, num);
    MPI_Type_commit(&gmc_row_type);

    // Initialize SIG matrices
    GMC_STEP("Init: SIG matrices");
    S0 = malloc(m * sizeof(sparse_matrix));
    ed = malloc(m * sizeof(matrix));
    sums = malloc(m * pattern_cnts[rank] * sizeof(double));
    __gmc_init_s0(X, m, num, S0, ed, sums);

    // U starts as average of SIG matrices
    GMC_STEP("Init: U");
    U = __gmc_init_u(S0, m, num);

    // Get matrix of eigenvectors (F), as well as eigenvalues
    GMC_STEP("Init: F");
    if(!rank) {
        F = new_matrix(num, num);
        F_old = new_matrix(num, num);
        evs = new_matrix(c + 1, NITER + 1);
    }

    MPI_Gatherv(U.data, counts[rank], gmc_row_type, F.data, counts, displs, gmc_row_type, 0, comm);
    F = update_f(F, evs.data, c, rank, blacs_row, blacs_col, blacs_height, blacs_width, blacs_ctx, comm, handle);

    // Initialize w to m uniform (All views start with the same weight)
    GMC_STEP("Init: w");
    double wI = 1.0 / m;
    w = new_matrix(m, 1);
    for(int v = 0; v < m; v++) w.data[v] = wI;

    // Main loop
    int it;
    for(it = 0; it < NITER; it++)
        if(__gmc_main_loop(it, S0, U, w, &F, &F_old, evs, m, c, num, ed, sums, &lambda))
            break;

    // Gather data to build final result
    matrix local_U = U;
    if(!rank) U = new_matrix(num, num);
    MPI_Gatherv(local_U.data, counts[rank], gmc_row_type, U.data, counts, displs, gmc_row_type, 0, comm);
    free_matrix(local_U);

    // Row type differs on sparse matrices
    MPI_Type_free(&gmc_row_type);
    gmc_row_type = gmc_contiguous_long(MPI_BYTE, (long long) (PN + 1) * sizeof(sprs_val));
    MPI_Type_commit(&gmc_row_type);

    // Collect SIG for each view
    for(int v = 0; v < m; v++) {
        sparse_matrix local_SIG = S0[v];
        if(!rank) S0[v] = new_sparse(PN + 1, num);
        MPI_Gatherv(local_SIG.data, counts[rank], gmc_row_type, S0[v].data, counts, displs, gmc_row_type, 0, comm);
        free_sparse(local_SIG);
    }

    // Cleanup, with workers
    MPI_Type_free(&gmc_row_type);
    for(int i = 0; i < m; i++) free_matrix(ed[i]);
    free(sums); free(ed); free_matrix(w);
    free(displs); free(counts); free(pattern_cnts);
    #ifdef ELPA_API_VER
        int error;
        elpa_deallocate(handle, &error);
        elpa_uninit(&error);
    #endif

    // Workers return before final clustering
    if(rank) {
        gmc_result null_result;
        return null_result;
    }

    // Adjacency matrix
    GMC_STEP("End: symU");
    bool * adj = malloc((long long) num * (long long) num * sizeof(bool));
    #pragma omp parallel for
    for(long long j = 0; j < num; j++)
        for(long long x = 0; x < j; x++)
            adj[j * num + x] = (U.data[j * num + x] != 0.0) || (U.data[x * num + j] != 0.0);

    // Final clustering. Find connected components on sU with Tarjan's algorithm
    GMC_STEP("End: final clustering");
    int * y = malloc((long long) num * sizeof(int));
    int cluster_num = connected_comp(adj, y, num);

    // Cleanup
    GMC_STEP("End: cleanup");
    free_matrix(F_old); free(adj);

    // Build output struct
    gmc_result result;
    result.U = U; result.S0 = S0; result.F = F; result.evs = evs; result.y = y; result.n = num; result.m = m;
    result.cluster_num = cluster_num; result.iterations = it + 1; result.lambda = lambda;
    return result;
}

void free_gmc_result(gmc_result r) {
    free_matrix(r.U);
    for(int i = 0; i < r.m; i++) free_sparse(r.S0[i]);
    free(r.S0);
    free_matrix(r.F);
    free_matrix(r.evs);
    free(r.y);
}
