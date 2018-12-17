#include "src/base.h"
#include "src/chebyshev.h"
#include "src/vector_operations.h"
#include "stdlib.h"

int main(int argc, char ** argv) {
    if (argc < 4) {
        printf("Usage:\n");
        printf("%s inputDataPath precision sParameter\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    const char * inputDataPath = argv[1];
    double precision = atof(argv[2]);
    int sParameter = atoi(argv[3]);
    int iterations = 0;

    struct Data data = loadDataFromFile(inputDataPath);
    int dimension = data.dimension;

    double * result = solveLinear(data, precision, sParameter, &iterations);

    printf("FINAL RESULT AFTER %i ITERATIONS IS:\n", iterations);
    printVector(result, dimension);

    deallocateData(data);
    free(result);

    return 0;
}
