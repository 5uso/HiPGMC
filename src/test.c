#include "gmc_scale.h"
#include "gmc_util.h"
#include "gmc_io.h"
#include "gmc.h"

#include <sys/time.h>
#include <stdio.h>
#include <mpi.h>

static dataset data[] = {
    {.path =    "../data/TwoMoon", .views = 2, .clusters =    2, .lambda = 1.0, .normalize = false},//0
    {.path =  "../data/ThreeRing", .views = 2, .clusters =    3, .lambda = 1.0, .normalize = false},//1
    {.path =        "../data/BBC", .views = 4, .clusters =    5, .lambda = 1.0, .normalize =  true},//2
    {.path =     "../data/Hdigit", .views = 2, .clusters =   10, .lambda = 1.0, .normalize =  true},//3
    {.path =  "../data/100leaves", .views = 3, .clusters =  100, .lambda = 1.0, .normalize =  true},//4
    {.path = "../data/Nuswide20k", .views = 5, .clusters =   81, .lambda = 1.0, .normalize = false},//5
    {.path =       "../data/Aloi", .views = 4, .clusters = 1000, .lambda = 1.0, .normalize =  true},//6
    {.path =    "../data/Airline", .views = 3, .clusters =    2, .lambda = 1.0, .normalize =  true},//7
    {.path =   "../data/Diabetes", .views = 1, .clusters =    2, .lambda = 1.0, .normalize =  true},//8
    {.path =      "../data/Covid", .views = 2, .clusters =    7, .lambda = 1.0, .normalize =  true},//9
    {.path =   "../data/3sources", .views = 3, .clusters =    6, .lambda = 1.0, .normalize =  true},//10
    {.path =   "../data/BBCSport", .views = 2, .clusters =    5, .lambda = 1.0, .normalize =  true},//11
    {.path =      "../data/Mfeat", .views = 6, .clusters =   10, .lambda = 1.0, .normalize =  true},//12
    {.path =        "../data/NGs", .views = 3, .clusters =    5, .lambda = 1.0, .normalize =  true},//13
    {.path =      "../data/WebKB", .views = 3, .clusters =    4, .lambda = 1.0, .normalize =  true},//14
};

int main(int argc, char *argv[]) {
    // MPI Initialization
    int rank, thread_support, numprocs;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &thread_support);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(thread_support != MPI_THREAD_MULTIPLE) {
        printf("MPI concurrent messaging not supported. Exiting\n");
        exit(1);
    }

    // BLACS Initialization
    int context, blacs_rows, blacs_cols;
    grid_dims(numprocs, &blacs_rows, &blacs_cols);
    if(!rank) printf("Process grid: %d rows, %d cols\n", blacs_rows, blacs_cols);
    Cblacs_get(0, 0, &context);
    Cblacs_gridinit(&context, "C", blacs_rows, blacs_cols);

    if(!rank) {
        dataset d = data[0];
        printf("Loading dataset from '%s'\n", d.path);
        matrix * X = read_dataset(d.path);
        //dataset d = {.path = "", .views = 3, .clusters = 4, .lambda = 1.0, .normalize = false};
        //matrix * X = generate_data(40000, 3, 100, 4, -1.0);
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
    
    Cblacs_gridexit(context);
    MPI_Finalize();
    return 0;
}
