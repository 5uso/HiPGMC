#include "gmc_sum.h"

double block_sum(double * arr, int n) {
    if(n == 0) return 0.0;

    double sum0 = 0.0, sum1 = 0.0;

    for(int i = 0; i < n - 1; i++)
        sum0 += arr[i], sum1 += arr[++i];

    return (sum0 + sum1) + arr[n - 1] * (n % 2);
}

double block_sum_col(double * arr, int n, int width) {
    if(n == 0) return 0.0;

    double sum0 = 0.0, sum1 = 0.0;
    long long w = width;

    for(int i = 0; i < n - 1; i++)
        sum0 += arr[i * w], sum1 += arr[++i * w];

    return (sum0 + sum1) + arr[(n - 1) * w] * (n % 2);
}

double block_sum_ptr(double ** arr, int n, long long offset) {
    if(n == 0) return 0.0;

    double sum0 = 0.0, sum1 = 0.0;

    for(int i = 0; i < n - 1; i++)
        sum0 += *(offset + arr[i]), sum1 += *(offset + arr[++i]);

    return (sum0 + sum1) + *(offset + arr[n - 1]) * (n % 2);
}

double block_sum_col_sqr(double * arr, int n, int width) {
    if(n == 0) return 0.0;

    double sum0 = 0.0, sum1 = 0.0;
    long long w = width;

    for(int i = 0; i < n - 1; i++) {
        sum0 += arr[i * w] * arr[i * w];
        i++;
        sum1 += arr[i * w] * arr[i * w];
    }

    return (sum0 + sum1) + arr[(n - 1) * w] * arr[(n - 1) * w] * (n % 2);
}
