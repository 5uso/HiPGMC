#ifndef GMCINPUT
#define GMCINPUT

#include "matrix.h"
#include <dirent.h>

matrix read_matrix(const char * path);
void dump_matrix(matrix m, const char * path);
matrix * read_dataset(const char * path);

#endif
