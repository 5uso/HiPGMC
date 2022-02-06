#include "heap.h"

void sift_down(heap h, uint i) {
    uint l = 2 * i;
    uint r = l + 1;

    uint max = i;
    if(l < h.size && *h.data[l] > *h.data[max]) max = l;
    if(r < h.size && *h.data[r] > *h.data[max]) max = r;

    if(max != i) {
        double * d = h.data[max];
        h.data[max] = h.data[i];
        h.data[i] = d;

        sift_down(h, max);
    }
}

heap newHeap(double * data, uint size) {
    heap h;
    h.data = malloc(sizeof(double*) * size);
    h.size = size;
    h.min  = data;

    for(int i = 0; i < size; i++) {
        h.data[i] = data++;
        if(*h.data[i] < *h.min) h.min = h.data[i];
    }
    for(int i = size / 2; i >= 0; i--) sift_down(h, i);

    return h;
}

void freeHeap(heap h) {
    free(h.data);
}

double heapMax(heap h) {
    return *h.data[0];
}

double heapMin(heap h) {
    return *h.min;
}

void replace(heap * h, double * val) {
    *h->data[0] = 0.0;
    h->data[0] = val;
    sift_down(*h, 0);
    if(*val < *h->min) h->min = val;
}