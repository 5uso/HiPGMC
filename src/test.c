#include "funs.h"
#include "heap.h"
#include "gmc.h"
#include "input.h"

#include <stdio.h>
#include <sys/time.h>

typedef struct dataset {
    const char * path;
    uint views;
    uint clusters;
    double lambda;
    bool normalize;
} dataset;

static dataset data[] = {
    {.path =   "../data/TwoMoon", .views = 2, .clusters =  2, .lambda = 1.0d, .normalize = false},
    {.path = "../data/ThreeRing", .views = 2, .clusters =  3, .lambda = 1.0d, .normalize = false},
    {.path =       "../data/BBC", .views = 4, .clusters =  5, .lambda = 1.0d, .normalize = true},
    {.path =    "../data/Hdigit", .views = 2, .clusters = 10, .lambda = 1.0d, .normalize = false},
};

void what_eig() {
    matrix F = read_matrix("../before_eig_o3");
    /*for(int y = 0; y < F.h; y++) {//L
        for(int x = 0; x < y - 2; x++)
            F.data[y * F.w + x] = 0.0d;
    }*/

    for(int y = 0; y < F.h; y++) {//U
        for(int x = y + 1; x < F.w; x++)
            F.data[y * F.w + x] = 0.0d;
    }

    int E = LAPACK_COL_MAJOR;
    int found_eigenvalue_n;
    int * ifail = malloc(sizeof(int) * 6);
    double * ev = malloc(sizeof(double) * 6);
    double * eigenvectors = malloc(sizeof(double) * F.w * 6);
    LAPACKE_dsyevx(E, 'V', 'I', 'U', F.w, F.data, F.w, -1.0d, -1.0d, 1, 6, 0.0d, &found_eigenvalue_n, ev, eigenvectors, F.w, ifail);
    F.h = 6;
    memcpy(F.data, eigenvectors, sizeof(double) * F.w * F.h);

    for(int i = 0; i < 6; i++) printf("%lf ", ev[i]);
    printf("\n");
}

void what_mul() {
    matrix m = read_matrix("../after_norm");
    matrix mt = new_matrix(m.w, m.w);
    //cblas_dsyrk(CblasRowMajor, CblasUpper, CblasTrans, m.w, m.h, 1.0d, m.data, m.w, 0.0d, mt.data, m.w);
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, m.w, m.w, m.h, 1.0d, m.data, m.w, m.data, m.w, 0.0d, mt.data, m.w);
    dump_matrix(mt, "mul_test2");
}

int main(int argc, char *argv[]) {
    /*what_eig();
    return 0;//*/
    /*what_mul();
    return 0;//*/

    dataset d = data[2];
    printf("Loading dataset from '%s'\n", d.path);
    matrix * X = read_dataset(d.path);
    if(!X) {
        perror("Couldn't load dataset");
        exit(1);
    }
    printf("Dataset loaded! Press enter to start\n");
    getchar();
    printf("Running GMC...\n");

    struct timeval begin, end;
    gettimeofday(&begin, 0);
    gmc_result r = gmc(X, d.views, d.clusters, d.lambda, d.normalize);
    gettimeofday(&end, 0);

    if(r.cluster_num != d.clusters) printf("Couldn't find requested cluster number (%d). Got %d clusters\n", d.clusters, r.cluster_num);
    printf("Iteration %d: λ=%lf\n", r.iterations, r.lambda);
    printf("Time measured: %.4f seconds.\n", end.tv_sec - begin.tv_sec + (end.tv_usec - begin.tv_usec) * 1e-6);
    
    for(int i = 0; i < r.n; i++) printf("%d ", r.y[i]);
    printf("\n");
    //print(r.U);

    free_gmc_result(r);

    return 0;
}
