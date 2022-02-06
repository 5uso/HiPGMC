#ifndef GMCHEAP
#define GMCHEAP

#include <stdlib.h>
#include <string.h>

typedef struct heap {
    double ** data;
    uint      size;
    double *   min;
} heap;

heap newHeap(double * data, uint size);
void freeHeap(heap h);
double heapMax(heap h);
double heapMin(heap h);
void replace(heap * h, double * val);

#endif