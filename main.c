#include "src/base.h"
#include "src/chebyshev.h"
#include "src/vector_operations.h"
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

    for (int sParameter = initSParameter; sParameter < initSParameter + calculationsNum; ++sParameter) {
        int iterations = 0;

        double *result = solveLinear(data, precision, sParameter, &iterations);

        free(result);
    }

    deallocateData(data);

    return 0;
}
