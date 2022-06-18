#include "gmc.h"

gmc_result gmc(matrix * X, uint m, uint c, double lambda, bool normalize) {
    uint num = X[0].w;

    //Normalize data
    if(normalize) {
        for(int i = 0; i < m; i++) {
            for(int y = 0; y < num; y++) {
                int w = X[i].w;
                double mean = 0.0d;
                for(int x = 0; x < w; x++) mean += X[i].data[y * w + x];
                mean /= (double) w;

                double std = 0.0d;
                for(int x = 0; x < w; x++) {
                    double dev = X[i].data[y * w + x] - mean;
                    std += dev * dev;
                }
                std /= (double) w;
                std = sqrt(std);

                for(int x = 0; x < w; x++) X[i].data[y * w + x] = (X[i].data[y * w + x] - mean) / (std + EPS);
            }
        }
    }

    //Initialize SIG matrices
    matrix * S0 = malloc(m * sizeof(matrix));
    for(int i = 0; i < m; i++) S0[i] = initSIG(X[i], PN);

    //U starts as average of SIG matrices
    matrix U = newMatrix(num, num);
    memcpy_s(U.data, num * num, S0[0].data, num * num);
    for(int i = 1; i < m; i++) {
        for(int y = 0; y < num; y++) {
            for(int x = 0; x < num; x++) U.data[y * num + x] += S0[i].data[y * num + x];
        }
    }
    for(int y = 0; y < num; y++) {
        for(int x = 0; x < num; x++) U.data[y * num + x] /= (double) num;
    }

    //Divide each row of U by its own sum
    for(int y = 0; y < num; y++) {
        double sum = 0.0d;
        for(int x = 0; x < num; x++) sum += U.data[y * num + x];
        for(int x = 0; x < num; x++) U.data[y * num + x] /= sum;
    }

    //Get matrix of eigenvectors (F), as well as eigenvalues
    matrix F = newMatrix(num, num);
    matrix F_old = newMatrix(num, num);
    matrix evs = newMatrix(num, NITER + 1);
    updateF(F, U, evs.data, c);

    // Initialize w to m uniform (All views start with the same weight)
    double wI = 1.0 / m;
    matrix w = newMatrix(m, 1);
    for(int i = 0; i < m; i++) w.data[i] = wI;

    // Used when calculating S0
    matrix * ed = malloc(m * sizeof(matrix));
    matrix * idxx = malloc(m * sizeof(matrix));
    for(int i = 0; i < m; i++) {
        ed[i] = sqrDist(X[i]);
        //TODO: Store sort into idxx (heap, since the loop uses lowest values?)
    }
    //After this is done we no longer need X, maybe free it.

    // Main loop
    for(int it = 0; it < NITER; it++) {
        //TODO: Update S0
        for(int i = 0; i < m; i++) {
            //S0 gets set to all zeros
            for(int y = 0; y < num; y++) {
                
            }
        }

        //Update w
        matrix US = newMatrix(num, num);
        for(int i = 0; i < m; i++) {
            memcpy_s(US.data, num * num, U.data, num * num);
            for(int y = 0; y < num; y++) {
                for(int x = 0; x < num; x++) US.data[y * num + x] -= S0[i].data[y * num + x];
            }

            double distUS = LAPACKE_dlange('F', num, num, US.data, num, NULL);
            w.data[i] = 0.5d / (distUS + EPS);
        }
        freeMatrix(US);

        //Update U
        matrix dist = sqrDist(F); //TODO: This actually needs F to be transposed
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
            matrix q = newMatrix(qw, m);
            for(int x = 0, i = 0; x < num; x++) {
                if(idx[x]) q.data[i++] = lambda * (double) abs(y - x) / (double) m * -0.5d;
            }
            for(int v = m - 1; v >= 0; v--) {
                for(int x = 0, i = 0; x < num; x++) {
                    if(idx[x]) q.data[v * qw + i] = q.data[i++] / w.data[v] + S0[v].data[y * num + x];
                }
            }
            
            q = updateU(q);
            for(int x = 0, i = 0; x < num; x++) {
                if(idx[x]) U.data[y * num + x] = q.data[i++];
                else U.data[y * num + x] = 0.0d;
            }
            freeMatrix(q);
        }
        free(idx);

        //Update matrix of eigenvectors (F), as well as eigenvalues
        matrix temp = F_old;
        F_old = F;
        F = temp;
        double * ev = evs.data + num * it;
        updateF(F, U, ev, c);

        //Update lambda
        double fn = 0.0d;
        for(int i = 0; i < c; i++) fn += ev.data[i];
        if(fn > ZR) {
            lambda *= 2.0d;
        } else if(fn + ev.data[c] < ZR) {
            lambda /= 2.0d;
            temp = F_old;
            F_old = F;
            F = temp;
        } else {
            printf("Iteration %d: λ=%lf\n", it, lambda);
            evs.h = it + 2;
            break;
        }
    }

    //TODO: Final clustering. Find connected components on sU with Tarjan's algorithm
    gmc_result result;
    result.y = NULL; result.U = U; result.S0 = S0; result.F = F; result.evs = evs;
    return result;
}