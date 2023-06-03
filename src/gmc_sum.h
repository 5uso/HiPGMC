#ifndef GMCSUM
#define GMCSUM

double block_sum(double * arr, int n);
double block_sum_col(double * arr, int n, int width);
double block_sum_ptr(double ** arr, int n, long long offset);
double block_sum_col_sqr(double * arr, int n, int width);

#endif
