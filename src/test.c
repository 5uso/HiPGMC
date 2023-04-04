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
    {.path =       "../data/BBC", .views = 4, .clusters =  5, .lambda = 1.0d, .normalize = false},
    {.path =    "../data/Hdigit", .views = 2, .clusters = 10, .lambda = 1.0d, .normalize = false},
};

int main(int argc, char *argv[]) {
    dataset d = data[0];
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
    
    //for(int i = 0; i < r.n; i++) printf("%d ", r.y[i]);
    //printf("\n");
    //print(r.U);

    free_gmc_result(r);

    return 0;
}
