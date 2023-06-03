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

GMC_INTERNAL matrix __gmc_init_u(matrix * S0, uint m, uint num) {
    matrix U = new_matrix(num, num);

    // U starts as average of SIG matrices
    memcpy(U.data, S0[0].data, num * num * sizeof(double));

    for(int v = 1; v < m; v++)
        for(int y = 0; y < num; y++)
            for(int x = 0; x < num; x++)
                U.data[y * num + x] += S0[v].data[y * num + x];

    for(int y = 0; y < num; y++)
        for(int x = 0; x < num; x++)
            U.data[y * num + x] /= m;

    // Divide each row of U by its own sum
    for(int y = 0; y < num; y++) {
        double sum = block_sum(U.data + y * num, num);
        for(int x = 0; x < num; x++) U.data[y * num + x] /= sum;
    }

    return U;
}

GMC_INTERNAL void __gmc_update_s0(matrix * S0, matrix U, matrix w, uint m, uint num, matrix * ed, heap * idxx, double * sums) {
    for(int v = 0; v < m; v++) {
        memset(S0[v].data, 0x00, num * num * sizeof(double));
        long long offsetU = U.data - ed[v].data;
        long long offsetS = S0[v].data - ed[v].data;
        double weight = w.data[v] * 2.0;

        for(int y = 0; y < num; y++) {
            heap h = idxx[v * num + y];
            double max = heap_max(h);
            double maxU = *(offsetU + h.data[0]);

            double sumU = block_sum_ptr(h.data + 1, PN, offsetU);

            double numerator = max - weight * maxU;
            double denominator = PN * max - sums[v * num + y] + weight * (sumU - PN * maxU) + EPS;

            for(int x = 0; x < PN + 1; x++) {
                double r = (numerator - *h.data[x] + weight * *(offsetU + h.data[x])) / denominator;
                *(offsetS + h.data[x]) = r * (r > 0.0);
            }
        }
    }
}

GMC_INTERNAL void __gmc_update_w(matrix * S0, matrix U, matrix w, uint m, uint num) {
    matrix US = new_matrix(num, num);

    for(int v = 0; v < m; v++) {
        memcpy(US.data, U.data, num * num * sizeof(double));
        for(int y = 0; y < num; y++)
            for(int x = 0; x < num; x++)
                US.data[y * num + x] -= S0[v].data[y * num + x];

        double distUS = LAPACKE_dlange(LAPACK_COL_MAJOR, 'F', num, num, US.data, num);
        w.data[v] = 0.5 / (distUS + EPS);
    }

    free_matrix(US);
}

GMC_INTERNAL void __gmc_update_u(matrix * S0, matrix U, matrix w, matrix * F, uint m, uint num, double * lambda) {
    matrix dist = sqr_dist(*F); // F is transposed, since LAPACK returns it in column major
    bool * idx = malloc(num * sizeof(bool));

    for(int y = 0; y < num; y++) {
        int qw = 0;

        #ifdef IS_LOCAL
            for(int x = 0; x < num; x++) {
                idx[x] = (S0[0].data[y * num + x] > 0);
                qw += idx[x];
            }

            for(int v = 1; v < m; v++) {
                for(int x = 0; x < num; x++) {
                    qw -= idx[x];
                    idx[x] |= (S0[v].data[y * num + x] > 0);
                    qw += idx[x];
                }
            }
        #else
            memset(idx, 0x01, num * sizeof(bool));
            qw = num;
        #endif

        matrix q = new_matrix(qw, m);
        for(int x = 0, i = 0; x < num; x++) {
            if(idx[x]) {
                q.data[i] = *lambda * dist.data[y * num + x] / (double) m * -0.5;
                i++;
            } 
        }

        for(int v = m - 1; v >= 0; v--) {
            for(int x = 0, i = 0; x < num; x++) {
                if(idx[x]) {
                    q.data[v * qw + i] = q.data[i] / w.data[v] + S0[v].data[y * num + x];
                    i++;
                }
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

GMC_INTERNAL bool __gmc_main_loop(int it, matrix * S0, matrix U, matrix w, matrix * F, matrix * F_old, matrix evs, uint m, uint c, uint num,
                                  matrix * ed, heap * idxx, double * sums, double * lambda) {
    // Update S0
    GMC_STEP(printf("Iteration %d, update S0\n", it));
    __gmc_update_s0(S0, U, w, m, num, ed, idxx, sums);

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
    matrix * S0 = malloc(m * sizeof(matrix));
    for(int v = 0; v < m; v++) S0[v] = init_sig(X[v], PN);

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

    // Used when calculating S0
    GMC_STEP(printf("Init, S0 sort\n"));
    matrix * ed = malloc(m * sizeof(matrix));
    heap * idxx = malloc(m * num * sizeof(heap));
    double * sums = malloc(m * num * sizeof(double));
    for(int v = 0; v < m; v++) {
        ed[v] = sqr_dist(X[v]);
        for(int y = 0; y < num; y++) {
            ed[v].data[y * num + y] = INFINITY;
            
            heap h = new_heap(ed[v].data + y * num, PN + 1);
            for(int x = PN + 1; x < num; x++)
                if(ed[v].data[y * num + x] < heap_max(h))
                    replace(&h, ed[v].data + y * num + x);
                
            idxx[v * num + y] = h;
            sums[v * num + y] = block_sum_ptr(h.data + 1, PN, 0);
        }
    }

    // Main loop
    int it;
    for(it = 0; it < NITER; it++)
        if(__gmc_main_loop(it, S0, U, w, &F, &F_old, evs, m, c, num, ed, idxx, sums, &lambda))
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
    for(int i = 0; i < m * num; i++) free_heap(idxx[i]);
    free_matrix(F_old); free_matrix(w);
    free(sums); free(idxx); free(ed); free(adj);

    // Build output struct
    gmc_result result;
    result.U = U; result.S0 = S0; result.F = F; result.evs = evs; result.y = y; result.n = num; result.m = m;
    result.cluster_num = cluster_num; result.iterations = it + 1; result.lambda = lambda;
    return result;
}

void free_gmc_result(gmc_result r) {
    free_matrix(r.U);
    for(int i = 0; i < r.m; i++) free_matrix(r.S0[i]);
    free(r.S0);
    free_matrix(r.F);
    free_matrix(r.evs);
    free(r.y);
}
