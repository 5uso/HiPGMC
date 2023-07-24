#ifndef GMCTEST
#define GMCTEST

#include "gmc_matrix.h"

#include <stdbool.h>

typedef struct dataset {
    const char * path;
    int views;
    int clusters;
    double lambda;
    bool normalize;
} dataset;

matrix * generate_data(int samples, int views, int features, int clusters, double difficulty);

#endif
