#include "gmc.h"

gmc_result gmc(matrix * X, uint m, uint c, double lambda, bool normalize) {
    uint num = X[0].w;

    //Normalize data
    if(normalize) {
        for(int i = 0; i < m; i++) {
            for(int y = 0; y < num; y++) {
                double mean = 0.0d;
                for(int x = 0; x < num; x++) mean += U.data[y][x];
                mean /= (double) num;

                double std = 0.0d;
                for(int x = 0; x < num; x++) {
                    double dev = U.data[y][x] - mean;
                    std += dev * dev;
                }
                std /= (double) num;
                std = sqrt(std);

                for(int x = 0; x < num; x++) U.data[y][x] = (U.data[y][x] - mean) / (std + EPS);
            }
        }
    }

    //Initialize SIG matrices
    matrix * S0 = malloc(m * sizeof(matrix *));
    for(int i = 0; i < m; i++) S0[i] = initSIG(X[i], PN);

    //U starts as average of SIG matrices
    matrix U = newMatrix(num, num);
    memcpy_s(U.data, num * num, S0[0].data, num * num);
    for(int i = 1; i < m; i++) {
        for(int y = 0; y < num; y++) {
            for(int x = 0; x < num; x++) U.data[y][x] += S0[i].data[y][x];
        }
    }
    for(int y = 0; y < num; y++) {
        for(int x = 0; x < num; x++) U.data[y][x] /= (double) num;
    }

    //Divide each row of U by its own sum
    for(int y = 0; y < num; y++) {
        double sum = 0.0d;
        for(int x = 0; x < num; x++) sum += U.data[y][x];
        for(int x = 0; x < num; x++) U.data[y][x] /= sum;
    }

    //TODO: The thing with symmetric U into eigenvectors, gets F
    matrix F;


    double wI = 1.0 / m;
    matrix w = newMatrix(m, 1);
    for(int i = 0; i < m; i++) w.data[i] = wI;

    matrix * ed = malloc(m * sizeof(matrix *));
    matrix * idxx = malloc(m * sizeof(matrix *));
    for(int i = 0; i < m; i++) {
        ed[i] = sqrDist(X[i]);
        //TODO: Store sort into idxx (heap, since the loop uses lowest values?)
    }

    for(int it = 0; it < NITER; it++) {
        //TODO: Update S0
        for(int i = 0; i < m; i++) {
            //S0 gets set to all zeros
            for(int y = 0; y < num; y++) {
                
            }
        }

        //Update w
        matrix sU = newMatrix(num, num);
        for(int i = 0; i < m; i++) {
            memcpy_s(sU.data, num * num, U.data, num * num);
            for(int y = 0; y < num; y++) {
                for(int x = 0; x < num; x++) sU.data[y][x] -= S0[i].data[y][x];
            }

            double distUS = LAPACKE_dlange('F', num, num, sU.data, num, NULL);
            w.data[i] = 0.5d / (distUS + EPS);
        }
        freeMatrix(sU);

        //TODO: Update U
        matrix dist = sqrDist(F); //This actually needs F to be transposed
        //U gets set to all zeros
        for(int y = 0; y < num; y++) {

        }

        //TODO: The thing with symmetric U into eigenvectors, gets F

        //TODO: Lambda stuff
    }

    //TODO: Final clustering
    gmc_result result;
    result.y = NULL; result.U = U; result.S0 = S0; result.F = NULL; result.evs = NULL;
    return result;
}