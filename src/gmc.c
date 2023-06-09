#include "gmc.h"

GMC_INTERNAL void __gmc_normalize(matrix * X, uint m, uint num) {
    for(int v = 0; v < m; v++) {
        int h = X[v].h, w = X[v].w;
        for(int x = 0; x < w; x++) {
            double mean = 0.0;
            for(int y = 0; y < h; y++) mean += X[v].data[y * w + x];
            mean /= h;

            double std = 0.0;
            for(int y = 0; y < h; y++) {
                double dev = X[v].data[y * w + x] - mean;
                std += dev * dev;
            }
            std /= h - 1;
            std = sqrt(std);
            if(std == 0) std = EPS;

            for(int y = 0; y < h; y++) X[v].data[y * w + x] = (X[v].data[y * w + x] - mean) / std;
        }
    }
}

GMC_INTERNAL void __gmc_init_s0(matrix * X, uint m, uint num, sparse_matrix * S0, matrix * ed, double * sums) {
    for(int v = 0; v < m; v++) {
        matrix ted = sqr_dist(X[v]);
        S0[v] = new_sparse(PN + 1, ted.h);
        ed[v] = new_matrix(PN + 1, ted.h);

        for(int y = 0; y < num; y++) {
            ted.data[y * num + y] = INFINITY;
            
            heap h = new_heap(ted.data + y * num, PN + 1);
            for(int x = PN + 1; x < num; x++)
                if(ted.data[y * num + x] < heap_max(h))
                    replace(&h, ted.data + y * num + x);
                
            sums[v * num + y] = block_sum_ptr(h.data + 1, PN, 0);
            double denominator = *h.data[0] * PN - sums[v * num + y] + EPS;

            for(int i = 0; i < PN + 1; i++) {
                sprs_val val = {
                    .i = h.data[i] - (ted.data + y * num),
                    .value = ((*h.data[0] - *h.data[i]) / denominator) * (i > 0),
                };
                S0[v].data[y * (PN + 1) + i] = val;
                ed[v].data[y * (PN + 1) + i] = *h.data[i];
            }

            free_heap(h);
        }

        free_matrix(ted);
    }
}

GMC_INTERNAL matrix __gmc_init_u(sparse_matrix * S0, uint m, uint num) {
    matrix U = new_matrix(num, num);

    // U starts as average of SIG matrices
    memset(U.data, 0x00, num * num * sizeof(double));

    for(int y = 0; y < num; y++) {
        double sum = 0.0;
        for(int i = 0; i < PN + 1; i++)
            for(int v = 0; v < m; v++) {
                sprs_val val = S0[v].data[y * (PN + 1) + i];
                double t = val.value / m;
                U.data[y * num + val.i] += t;
                sum += t;
            }
        for(int x = 0; x < num; x++) U.data[y * num + x] /= sum;
    }

    return U;
}

GMC_INTERNAL void __gmc_update_s0(sparse_matrix * S0, matrix U, matrix w, uint m, uint num, matrix * ed, double * sums) {
    for(int v = 0; v < m; v++) {
        double weight = w.data[v] * 2.0;

        for(int y = 0; y < num; y++) {
            double max = ed[v].data[(PN + 1) * y];
            double maxU = U.data[y * num + S0[v].data[(PN + 1) * y].i];

            double sumU = 0.0;
            for(int i = 1; i < PN + 1; i++) {
                int x = S0[v].data[y * (PN + 1) + i].i;
                sumU += U.data[y * num + x];
            }

            double numerator = max - weight * maxU;
            double denominator = PN * max - sums[v * num + y] + weight * (sumU - PN * maxU) + EPS;

            for(int i = 0; i < PN + 1; i++) {
                int x = S0[v].data[y * (PN + 1) + i].i;
                double r = (numerator - ed[v].data[(PN + 1) * y + i] + weight * U.data[y * num + x]) / denominator;
                S0[v].data[y * (PN + 1) + i].value = r * (r > 0.0);
            }
        }
    }
}

GMC_INTERNAL void __gmc_update_w(sparse_matrix * S0, matrix U, matrix w, uint m, uint num) {
    matrix US = new_matrix(num, num);

    for(int v = 0; v < m; v++) {
        memcpy(US.data, U.data, num * num * sizeof(double));
        for(int y = 0; y < num; y++)
            for(int i = 0; i < PN + 1; i++) {
                sprs_val val = S0[v].data[y * (PN + 1) + i];
                int x = val.i;
                US.data[y * num + x] -= val.value;
            }

        double distUS = LAPACKE_dlange(LAPACK_COL_MAJOR, 'F', num, num, US.data, num);
        w.data[v] = 0.5 / (distUS + EPS);
    }

    free_matrix(US);
}

GMC_INTERNAL void __gmc_update_u(sparse_matrix * S0, matrix U, matrix w, matrix * F, uint m, uint num, double * lambda) {
    matrix dist = sqr_dist(*F); // F is transposed, since LAPACK returns it in column major
    int * idx = malloc(num * sizeof(int));

    for(int y = 0; y < num; y++) {
        int qw = 0;

        #ifdef IS_LOCAL
            memset(idx, 0x00, num * sizeof(int));
            for(int v = 0; v < m; v++) {
                for(int i = 0; i < PN + 1; i++) {
                    sprs_val val = S0[v].data[y * (PN + 1) + i];
                    int x = val.i;
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
        for(int x = 0, i = 0; x < num; x++) {
            if(idx[x]) {
                q.data[i] = *lambda * dist.data[y * num + x] / (double) m * -0.5;
                i++;
                idx[x] = i;
            } 
        }

        for(int v = m - 1; v >= 0; v--)
            for(int i = 0; i < qw; i++)
                q.data[v * qw + i] = q.data[i] / w.data[v];

        for(int v = m - 1; v >= 0; v--) {
            for(int i = 0; i < PN + 1; i++) {
                sprs_val val = S0[v].data[y * (PN + 1) + i];
                int x = val.i;
                if(idx[x]) q.data[v * qw + idx[x] - 1] += val.value;
            }
        }

        q = update_u(q);
        for(int x = 0, i = 0; x < num; x++) {
            if(idx[x]) {
                U.data[y * num + x] = q.data[i];
                i++;
            } else U.data[y * num + x] = 0.0;
        }

        free_matrix(q);
    }

    free(idx);
    free_matrix(dist);
}

GMC_INTERNAL bool __gmc_main_loop(int it, sparse_matrix * S0, matrix U, matrix w, matrix * F, matrix * F_old, matrix evs, uint m, uint c, uint num,
                                  matrix * ed, double * sums, double * lambda) {
    // Update S0
    GMC_STEP(printf("Iteration %d, update S0\n", it));
    __gmc_update_s0(S0, U, w, m, num, ed, sums);

    // Update w
    GMC_STEP(printf("Iteration %d, update w\n", it));
    __gmc_update_w(S0, U, w, m, num);

    // Update U
    GMC_STEP(printf("Iteration %d, update U\n", it));
    __gmc_update_u(S0, U, w, F, m, num, lambda);

    // Update matrix of eigenvectors (F), as well as eigenvalues
    GMC_STEP(printf("Iteration %d, update F\n", it));
    matrix temp = *F_old;
    *F_old = *F;
    *F = temp;
    double * ev = evs.data + (c + 1) * it;
    *F = update_f(*F, U, ev, c);

    // Update lambda
    GMC_STEP(printf("Iteration %d, update lambda\n", it));
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
        return true;
    }

    return false;
}

gmc_result gmc(matrix * X, uint m, uint c, double lambda, bool normalize) {
    uint num = X[0].w;

    // Normalize data
    GMC_STEP(printf("Init, normalize\n"));
    if(normalize) __gmc_normalize(X, m, num);

    // Initialize SIG matrices
    GMC_STEP(printf("Init, SIG matrices\n"));
    sparse_matrix * S0 = malloc(m * sizeof(sparse_matrix));
    matrix * ed = malloc(m * sizeof(matrix));
    double * sums = malloc(m * num * sizeof(double));
    __gmc_init_s0(X, m, num, S0, ed, sums);

    // U starts as average of SIG matrices
    GMC_STEP(printf("Init, U\n"));
    matrix U = __gmc_init_u(S0, m, num);

    // Get matrix of eigenvectors (F), as well as eigenvalues
    GMC_STEP(printf("Init, F\n"));
    matrix F = new_matrix(num, num);
    matrix F_old = new_matrix(num, num);
    matrix evs = new_matrix(c + 1, NITER + 1);
    F = update_f(F, U, evs.data, c);

    // Initialize w to m uniform (All views start with the same weight)
    GMC_STEP(printf("Init, w\n"));
    double wI = 1.0 / m;
    matrix w = new_matrix(m, 1);
    for(int v = 0; v < m; v++) w.data[v] = wI;

    // Main loop
    int it;
    for(it = 0; it < NITER; it++)
        if(__gmc_main_loop(it, S0, U, w, &F, &F_old, evs, m, c, num, ed, sums, &lambda))
            break;

    // Adjacency matrix
    GMC_STEP(printf("End, symU\n"));
    bool * adj = malloc(num * num * sizeof(bool));
    for(int j = 0; j < num; j++)
        for(int x = 0; x < j; x++)
            adj[j * num + x] = (U.data[j * num + x] != 0.0) || (U.data[x * num + j] != 0.0);

    // Final clustering. Find connected components on sU with Tarjan's algorithm
    GMC_STEP(printf("End, final clustering\n"));
    int * y = malloc(num * sizeof(int));
    int cluster_num = connected_comp(adj, y, num);

    // Cleanup
    GMC_STEP(printf("End, cleanup\n"));
    for(int i = 0; i < m; i++) free_matrix(ed[i]);
    free_matrix(F_old); free_matrix(w);
    free(sums); free(ed); free(adj);

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
