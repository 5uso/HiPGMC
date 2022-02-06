#include "funs.h"
#include "heap.h"

#include <stdio.h>
#include <sys/time.h>

#include <lapacke.h>

int main(int argc, char *argv[]) {
    struct timeval begin, end; matrix m;
    gettimeofday(&begin, 0);
    m = newMatrix(3, 3);
    set(m, 0, 0, 2.0d);
    set(m, 1, 0, 6.0d);
    set(m, 2, 0, 10.0d);
    set(m, 0, 1, 6.0d);
    set(m, 1, 1, 10.0d);
    set(m, 2, 1, 14.0d);
    set(m, 0, 2, 10.0d);
    set(m, 1, 2, 14.0d);
    set(m, 2, 2, 18.0d);
    
    matrix d = eig(m, 3);
    print(d);
    printf("\n");
    print(m);
    gettimeofday(&end, 0);
    printf("Time measured: %.3f seconds.\n", end.tv_sec - begin.tv_sec + (end.tv_usec - begin.tv_usec) * 1e-6);

    return 0;
}