#include "gmc_io.h"
#include "gmc.h"

#include <sys/time.h>
#include <stdio.h>

typedef struct dataset {
    const char * path;
    uint views;
    uint clusters;
    double lambda;
    bool normalize;
} dataset;

static dataset data[] = {
    {.path =   "../data/TwoMoon", .views = 2, .clusters =   2, .lambda = 1.0, .normalize = false},//0
    {.path = "../data/ThreeRing", .views = 2, .clusters =   3, .lambda = 1.0, .normalize = false},//1
    {.path =       "../data/BBC", .views = 4, .clusters =   5, .lambda = 1.0, .normalize =  true},//2
    {.path =    "../data/Hdigit", .views = 2, .clusters =  10, .lambda = 1.0, .normalize =  true},//3
    {.path = "../data/100leaves", .views = 3, .clusters = 100, .lambda = 1.0, .normalize =  true},//4
    {.path =  "../data/3sources", .views = 3, .clusters =   6, .lambda = 1.0, .normalize =  true},//5
    {.path =  "../data/BBCSport", .views = 2, .clusters =   5, .lambda = 1.0, .normalize =  true},//6
    {.path =     "../data/Mfeat", .views = 6, .clusters =  10, .lambda = 1.0, .normalize =  true},//7
    {.path =       "../data/NGs", .views = 3, .clusters =   5, .lambda = 1.0, .normalize =  true},//8
    {.path =     "../data/WebKB", .views = 3, .clusters =   4, .lambda = 1.0, .normalize =  true},//9
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
    
    for(int i = 0; i < r.n; i++) printf("%d ", r.y[i]);
    printf("\n");
    //print(r.U);

    free_gmc_result(r);
    for(int i = 0; i < d.views; i++) free_matrix(X[i]);
    free(X);

    return 0;
}
