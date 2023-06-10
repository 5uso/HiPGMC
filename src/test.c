#include "gmc_scalapack.h"
#include "gmc_io.h"
#include "gmc.h"

#include <sys/time.h>
#include <stdio.h>
#include <mpi.h>

typedef struct dataset {
    const char * path;
    int views;
    int clusters;
    double lambda;
    bool normalize;
} dataset;

static dataset data[] = {
    {.path =   "../data/TwoMoon", .views = 2, .clusters =   2, .lambda = 1.0, .normalize = false},
    {.path = "../data/ThreeRing", .views = 2, .clusters =   3, .lambda = 1.0, .normalize = false},
    {.path =       "../data/BBC", .views = 4, .clusters =   5, .lambda = 1.0, .normalize =  true},
    {.path =    "../data/Hdigit", .views = 2, .clusters =  10, .lambda = 1.0, .normalize =  true},
    {.path = "../data/100leaves", .views = 3, .clusters = 100, .lambda = 1.0, .normalize =  true},
};

int main(int argc, char *argv[]) {
    // MPI Initialization
    int rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // BLACS Initialization
    int context, blacs_rows = 1, blacs_cols = 1;
    Cblacs_get(0, 0, &context);
    Cblacs_gridinit(&context, "C", blacs_rows, blacs_cols);

    if(!rank) {
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
        gmc_result r = gmc(X, d.views, d.clusters, d.lambda, d.normalize, MPI_COMM_WORLD, context);
        gettimeofday(&end, 0);

        if(r.cluster_num != d.clusters) printf("Couldn't find requested cluster number (%d). Got %d clusters\n", d.clusters, r.cluster_num);
        printf("Iteration %d: λ=%lf\n", r.iterations, r.lambda);
        printf("Time measured: %.4f seconds.\n", end.tv_sec - begin.tv_sec + (end.tv_usec - begin.tv_usec) * 1e-6);
        
        for(int i = 0; i < r.n; i++) printf("%d ", r.y[i]);
        printf("\n");

        free_gmc_result(r);
        for(int i = 0; i < d.views; i++) free_matrix(X[i]);
        free(X);
    } else {
        // Only rank 0 takes care of IO, other processes receive info from it
        gmc(NULL, 0, 0, 0, 0, MPI_COMM_WORLD, context);
    } 
    
    MPI_Finalize();
    return 0;
}
