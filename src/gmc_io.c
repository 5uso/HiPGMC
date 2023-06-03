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

    uint w, h;
    fscanf(fd, "%u %u", &w, &h);
    matrix m = new_matrix(w, h);

    for(int y = 0; y < m.h; y++)
        for(int x = 0; x < m.w; x++)
            fscanf(fd, "%lg", &m.data[y * m.w + x]);

    return m;
}

void dump_matrix(matrix m, const char * path) {
    FILE * fd = fopen(path, "w");

    fprintf(fd, "%u %u\n", m.w, m.h);
    for(int y = 0; y < m.h; y++) {
        for(int x = 0; x < m.w - 1; x++)
            fprintf(fd, "%.17lg ", m.data[y * m.w + x]);
        fprintf(fd, "%.17lg\n", m.data[y * m.w + m.w - 1]);
    }
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

    return data;
}
