#include "gmc.h"

gmc_result gmc(matrix * X, uint m, uint c, double lambda, bool normalize) {
    uint num = X[0].w;

    //Normalize data
    if(normalize) {
        //TODO
    }

    //Initialize SIG matrices
    matrix * S0 = malloc(m * sizeof(matrix *));
    for(int i = 0; i < m; i++) S0[i] = initSIG(X[i], PN);

    matrix U = newMatrix(num, num);
    //TODO: U starts as average of SIG matrices
    //TODO: Divide each row of U by its own sum

    //TODO: The thing with symmetric U into eigenvectors, gets F

    double wI = 1.0 / m;
    matrix w = newMatrix(m, 1);
    for(int i = 0; i < m; i++) w.data[i] = wI;

    matrix * ed = malloc(m * sizeof(matrix *));
    matrix * idxx = malloc(m * sizeof(matrix *));
    for(int i = 0; i < m; i++) {
        ed[i] = sqrDist(X[i]);
        //TODO: Store sort into idxx
    }

    for(int it = 0; i < NITER; i++) {
        //TODO: Update S0

        //TODO: Update w

        //TODO: Update U

        //TODO: The thing with symmetric U into eigenvectors, gets F

        //TODO: Lambda stuff
    }

    //TODO: Final clustering
}