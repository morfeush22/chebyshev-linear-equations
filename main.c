#include "src/base.h"
#include "src/chebyshev.h"
#include "src/vector_operations.h"
#include "mpi.h"
#include "stdlib.h"
#include "time.h"

int main(int argc, char ** argv) {
    int rank, size;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc < 5) {
        printf("Usage:\n");
        printf("%s inputDataPath precision initSParameter calculationsNum\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    const char * inputDataPath = argv[1];
    double precision = atof(argv[2]);
    int initSParameter = atoi(argv[3]);
    int calculationsNum = atoi(argv[4]);

    struct Data data;
    if (rank == 0) {
        data = loadDataFromFile(inputDataPath);
    }
    else {
        data.matrix = NULL;
        data.bVector = NULL;
        data.dimension = 0;
    }

    int iterations;

    struct timespec start, stop;
    unsigned long deltaUS;

    for (int sParameter = initSParameter; sParameter < initSParameter + calculationsNum; ++sParameter) {
        if (rank == 0) {
            clock_gettime(CLOCK_MONOTONIC_RAW, &start);
        }

        double *result = solveLinear(data, precision, sParameter, &iterations, rank, size);

        if (rank == 0) {
            clock_gettime(CLOCK_MONOTONIC_RAW, &stop);

            deltaUS = (unsigned long)
                              (stop.tv_sec - start.tv_sec) * 1000000 + (stop.tv_nsec - start.tv_nsec) / 1000;

            printf("%i %lu\n", sParameter, deltaUS);
        }

        free(result);
    }

    deallocateData(data);

    MPI_Finalize();

    return 0;
}
