#include "gmc.h"

GMC_INTERNAL void __gmc_normalize(matrix * X, uint m, uint num) {
    for(int v = 0; v < m; v++) {
        for(int y = 0; y < num; y++) {
            int w = X[v].w;
            double mean = 0.0d;
            for(int x = 0; x < w; x++) mean += X[v].data[y * w + x];
            mean /= (double) w;

            double std = 0.0d;
            for(int x = 0; x < w; x++) {
                double dev = X[v].data[y * w + x] - mean;
                std += dev * dev;
            }
            std /= (double) w;
            std = sqrt(std);

            for(int x = 0; x < w; x++) X[v].data[y * w + x] = (X[v].data[y * w + x] - mean) / (std + EPS);
        }
    }
}

GMC_INTERNAL matrix __gmc_init_u(matrix * S0, uint m, uint num) {
    matrix U = new_matrix(num, num);

    // U starts as average of SIG matrices
    memcpy(U.data, S0[0].data, num * num * sizeof(double));

    for(int v = 1; v < m; v++) {
        for(int y = 0; y < num; y++) {
            for(int x = 0; x < num; x++) U.data[y * num + x] += S0[v].data[y * num + x];
        }
    }

    for(int y = 0; y < num; y++) {
        for(int x = 0; x < num; x++) U.data[y * num + x] /= (double) num;
    }

    // Divide each row of U by its own sum
    for(int y = 0; y < num; y++) {
        double sum = 0.0d;
        for(int x = 0; x < num; x++) sum += U.data[y * num + x];
        for(int x = 0; x < num; x++) U.data[y * num + x] /= sum;
    }

    return U;
}

GMC_INTERNAL void __gmc_update_s0(matrix * S0, matrix U, matrix w, uint m, uint num, matrix * ed, heap * idxx, double * sums) {
    for(int v = 0; v < m; v++) {
        memset(S0[v].data, 0x00, num * num * sizeof(double));
        for(int y = 0; y < num; y++) {
            heap h = idxx[v * num + y];
            long long offsetU = U.data - ed[v].data;
            long long offsetS = S0[v].data - ed[v].data;
            double weight = w.data[v] *  2.0d;

            double sumU = -*(offsetU + h.min);
            for(int x = 1; x < PN + 2; x++) sumU += *(offsetU + h.data[x]);

            double numerator = heap_max(h) - *(offsetU + h.data[0]) * weight;
            double denominator1 = PN * heap_max(h) - sums[v * num + y];
            double denominator2 = (sumU - *(offsetU + h.data[0]) * PN) * weight;

            for(int x = 0; x < PN + 2; x++) {
                if(h.data[x] == h.min) continue;
                double r = (numerator - *h.data[x] + weight * *(offsetU + h.data[x])) / (denominator1 + denominator2 + EPS);
                if(r > 0.0d) *(offsetS + h.data[x]) = r;
            }
        }
    }
}

GMC_INTERNAL void __gmc_update_w(matrix * S0, matrix U, matrix w, uint m, uint num) {
    matrix US = new_matrix(num, num); // TODO: We don't need to reallocate this every time

    for(int v = 0; v < m; v++) {
        memcpy(US.data, U.data, num * num * sizeof(double));
        for(int y = 0; y < num; y++) {
            for(int x = 0; x < num; x++) US.data[y * num + x] -= S0[v].data[y * num + x];
        }

        double distUS = LAPACKE_dlange(LAPACK_COL_MAJOR, 'F', num, num, US.data, num);
        w.data[v] = 0.5d / (distUS + EPS);
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
            memset(idx, 0x00, num * sizeof(bool));
            qw = num;
        #endif

        matrix q = new_matrix(qw, m);
        for(int x = 0, i = 0; x < num; x++) {
            if(idx[x]) {
                q.data[i] = *lambda * dist.data[y * num + x] / (double) m * -0.5d;
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
            } else U.data[y * num + x] = 0.0d;
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
    double * ev = evs.data + num * it;
    *F = update_f(*F, U, ev, c);

    // Update lambda
    GMC_STEP(printf("Iteration %d, update lambda\n", it));
    double fn = 0.0d;
    for(int i = 0; i < c; i++) fn += ev[i];
    if(fn > ZR) {
        *lambda *= 2.0d;
    } else if(fn + ev[c] < ZR) {
        *lambda /= 2.0d;
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
    matrix evs = new_matrix(num, NITER + 1);
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
        // TODO: Check -> Store sort into idxx (heap, since the loop uses lowest values?)
        for(int y = 0; y < num; y++) {
            heap h = new_heap(ed[v].data + y * num, PN + 2);
            for(int x = PN + 2; x < num; x++) {
                if(ed[v].data[y * num + x] < heap_max(h)) replace(&h, ed[v].data + y * num + x);
            }
            idxx[v * num + y] = h;
            sums[v * num + y] = -heap_min(h);
            for(int x = 1; x < PN + 2; x++) sums[v * num + y] += *h.data[x];
        }
    }

    // Main loop
    int it;
    for(it = 0; it < NITER; it++) {
       if(__gmc_main_loop(it, S0, U, w, &F, &F_old, evs, m, c, num, ed, idxx, sums, &lambda)) break;
    }

    // U symmetric
    GMC_STEP(printf("End, symU\n"));
    matrix sU = new_matrix(num, num);
    for(int y = 0; y < num; y++) {
        for(int x = y + 1; x < num; x++) sU.data[y * num + x] = (U.data[y * num + x] + U.data[x * num + y]) / 2.0d;
    }

    // Final clustering. Find connected components on sU with Tarjan's algorithm
    GMC_STEP(printf("End, final clustering\n"));
    int * y = malloc(sU.w * sizeof(int));
    int cluster_num = connected_comp(sU, y);

    // Cleanup
    GMC_STEP(printf("End, cleanup\n"));
    for(int i = 0; i < m; i++) free_matrix(ed[i]);
    for(int i = 0; i < m * num; i++) free_heap(idxx[i]);
    free_matrix(F_old); free_matrix(sU); free_matrix(w);
    free(sums); free(idxx); free(ed);

    // Build output struct
    gmc_result result;
    result.U = U; result.S0 = S0; result.F = F; result.evs = evs; result.y = y; result.n = sU.w; result.m = m;
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
