#ifndef GMCHEAP
#define GMCHEAP

#include <stdlib.h>
#include <string.h>

typedef struct heap {
    double ** data;
    uint      size;
    double *   min;
} heap;

heap new_heap(double * data, uint size);
void free_heap(heap h);
double heap_max(heap h);
double heap_min(heap h);
void replace(heap * h, double * val);
double * heap_pop(heap * h);

#endif
