#include "gmc_io.h"

int __count_dir_files(const char * path) {
    int cnt = 0;
    struct dirent * entry;
    DIR * d = opendir(path);

    if(!d) return 0;

    while((entry = readdir(d))) cnt += (entry->d_type == DT_REG);

    closedir(d);
    return cnt;
}

matrix read_matrix(const char * path) {
    FILE * fd = fopen(path, "r");

    int w, h;
    fscanf(fd, "%u %u", &w, &h);
    matrix m = new_matrix(w, h);

    for(long long y = 0; y < m.h; y++)
        for(long long x = 0; x < m.w; x++)
            fscanf(fd, "%lg", &m.data[y * m.w + x]);

    fclose(fd);
    return m;
}

void dump_matrix(matrix m, const char * path) {
    FILE * fd = fopen(path, "w");

    fprintf(fd, "%u %u\n", m.w, m.h);
    for(long long y = 0; y < m.h; y++) {
        for(long long x = 0; x < m.w - 1; x++)
            fprintf(fd, "%.17lg ", m.data[y * m.w + x]);
        fprintf(fd, "%.17lg\n", m.data[y * m.w + m.w - 1]);
    }

    fclose(fd);
}

void dump_sparse(sparse_matrix m, const char * path, int width) {
    matrix dense = new_matrix(width, m.h);
    memset(dense.data, 0x00, width * m.h * sizeof(double));

    for(long long y = 0; y < m.h; y++) 
        for(int i = 0; i < m.w; i++) {
            sprs_val val = m.data[y * m.w + i];
            long long x = val.i;
            dense.data[y * width + x] = val.value;
        }

    dump_matrix(dense, path);
    free_matrix(dense);
}

matrix * read_dataset(const char * path) {
    int views = __count_dir_files(path);
    if(!views) return NULL;

    matrix * data = malloc(sizeof(matrix) * views);

    struct dirent * entry;
    DIR * d = opendir(path);

    for(int i = 0; i < views; i++) {
        entry = readdir(d);
        if(entry->d_type != DT_REG) {
            i--;
            continue;
        }

        char filepath[1024];
        snprintf(filepath, 1024, "%s/%s", path, entry->d_name);
        filepath[1023] = '\0';

        data[i] = read_matrix(filepath);
    }

    closedir(d);
    return data;
}
