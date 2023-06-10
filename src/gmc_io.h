#ifndef GMCIO
#define GMCIO

#include "gmc_matrix.h"

#include <string.h>
#include <dirent.h>

matrix read_matrix(const char * path);
void dump_matrix(matrix m, const char * path);
void dump_sparse(sparse_matrix m, const char * path);
matrix * read_dataset(const char * path);

#endif
