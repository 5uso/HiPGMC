<h1 align=center>
  <img src="https://github.com/5uso/HiPGMC/assets/41950530/816b2656-ed13-4d98-93a9-f88b000a3ce6" alt="HiPGMC logo" width="600"/><br><br>
  Highly efficient, parallel Graph-based Multi-view Clustering
</h1>

<p align="justify">
  <b>Hi</b>gh <b>P</b>erformance <b>G</b>raph-based <b>M</b>ulti-view <b>C</b>lustering (HiPGMC) is an implementation of the novel <a href="https://github.com/cshaowang/gmc">GMC</a> clustering algorithm exploiting both shared and distributed memory parallelism.
  It's written in C with a focus in simplicity and a lean footprint.
</p>

## Dependencies

To build and run HiPGMC, the following prerequisites must be met:
- **GNU make:** Build management is handled via a very simple makefile.
- **C compiler:** Both GNU’s gcc and Intel’s icc have been confirmed to work.
- **MPI:** Both OpenMPI and Intel MPI have been tested. A threading level of `MPI_THREAD_MULTIPLE` is required.
- **OpenMP runtime:** Both GNU and Intel alternatives work.
- **BLAS:** OpenBLAS recommended. For performance reasons, ensure your BLAS implementation is threaded. Intel’s MKL is also confirmed to work.
- **ScaLAPACK:** The reference implementation was tested. BLACS should be either included within ScaLAPACK or provided separately. Intel’s MKL also works.
- **ELPA:** ELPA should be built with OpenMP support and linked to the same library implementations as HiPGMC. For increased performance, make sure SIMD kernels are enabled.

## Building from sources

After cloning the repository, navigate to its root folder and run `make` with default settings to build the shared library. The output `.so` and header files will be in the `out`
directory. To link to the library, these should be added to the library and include paths, respectively. `make test` outputs a statically linked test executable.

## Example application

```c
#include <gmc_util.h>
#include <sys/time.h>
#include <stdio.h>
#include <mpi.h>
#include <gmc.h>

int main(int argc, char *argv[]) {
    // MPI Initialization
    int rank, ts, numprocs;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &ts);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(ts != MPI_THREAD_MULTIPLE) {
        printf("MPI concurrent messaging not supported. Exiting\n");
        exit(1);
    }

    // BLACS Initialization
    int context, blacs_rows, blacs_cols;
    grid_dims(numprocs, &blacs_rows, &blacs_cols);
    Cblacs_get(0, 0, &context);
    Cblacs_gridinit(&context, "C", blacs_rows, blacs_cols);

    if(!rank) {
        printf("Process grid: %d rows, %d cols\n", blacs_rows, blacs_cols);

        // Generate a toy dataset for testing
        dataset d = { .path = "", .views = 3, .clusters = 4, .lambda = 1, .normalize = 0 };
        matrix * X = generate_data(40000, 3, 100, 4, 0.7);
        printf("Dataset generated! Press enter to start\n");
        getchar();
        printf("Running GMC...\n");

        // Run GMC
        struct timeval begin, end;
        gettimeofday(&begin, 0);
        gmc_result r = gmc(X, d.views, d.clusters, d.lambda,
        d.normalize, MPI_COMM_WORLD, context);
        gettimeofday(&end, 0);

        // Display results
        if(r.cluster_num != d.clusters) printf("Couldn't find %d clusters. Got %d instead\n", d.clusters, r.cluster_num);
        printf("Iteration %d: λ=%lf\n", r.iterations, r.lambda);
        printf("Time measured: %.4f seconds.\n", end.tv_sec - begin.tv_sec + (end.tv_usec - begin.tv_usec) * 1e-6);

        // Clustering vector
        for(int i = 0; i < r.n; i++) printf("%d ", r.y[i]);
        printf("\n");

        // Cleanup
        free_gmc_result(r);
        for(int i = 0; i < d.views; i++) free_matrix(X[i]);
        free(X);
    } else {
        // Only rank 0 takes care of IO
        gmc(NULL, 0, 0, 0, 0, MPI_COMM_WORLD, context);
    }

    Cblacs_gridexit(context);
    MPI_Finalize();
    return 0;
}
```
