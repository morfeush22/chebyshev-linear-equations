#include "src/base.h"
#include "src/chebyshev.h"
#include "src/vector_operations.h"
#include "sys/time.h"
#include "stdlib.h"

int main(int argc, char ** argv) {
    if (argc < 5) {
        printf("Usage:\n");
        printf("%s inputDataPath precision initSParameter calculationsNum\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    const char * inputDataPath = argv[1];
    double precision = atof(argv[2]);
    int initSParameter = atoi(argv[3]);
    int calculationsNum = atoi(argv[4]);
    
    struct Data data = loadDataFromFile(inputDataPath);

    int iterations;

    struct timespec start, stop;
    unsigned long deltaUS;

    for (int sParameter = initSParameter; sParameter < initSParameter + calculationsNum; ++sParameter) {
        clock_gettime(CLOCK_MONOTONIC_RAW, &start);
        double *result = solveLinear(data, precision, sParameter, &iterations);
        clock_gettime(CLOCK_MONOTONIC_RAW, &stop);

        deltaUS = (unsigned long)
                (stop.tv_sec - start.tv_sec) * 1000000 + (stop.tv_nsec - start.tv_nsec) / 1000;

        printf("%i %lu\n", sParameter, deltaUS);

        free(result);
    }

    deallocateData(data);

    return 0;
}
