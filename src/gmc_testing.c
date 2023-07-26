#include "gmc_testing.h"

static inline double rand_d(double min, double max) {
    return min + (rand() / (RAND_MAX / (max - min)));
}

// Generate random data matrix that results in a cyclic clustering
matrix * generate_data(int samples, int views, int features, int clusters, double difficulty) {
    srand(0);
    matrix * r = malloc(views * sizeof(matrix));
    for(int v = 0; v < views; v++) {
        r[v] = new_matrix(samples, features);
        int pos = 0;
        for(int y = 0; y < features; y++) {
            int part = y % clusters;
            for(int x = 0; x < samples; x++, pos++) {
                int cluster = x % clusters;
                r[v].data[pos] = (part == cluster) && (rand_d(0.0, 1.0) > difficulty) ? rand_d(0.05, 1.0) : 0.0;
            }
        }
    }

    return r;
}
