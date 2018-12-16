#include "src/base.h"
#include "src/chebyshev.h"
#include "src/vector_operations.h"
#include "stdlib.h"

int main() {
    int sParameter = 100;
    int iterations = 20;

    struct Data data = loadDataFromFile("tests/sources/simple_eq");

    double * result = solveLinear(data, sParameter, iterations);

    printf("FINAL RESULT IS:\n");
    printVector(result, data.dimension);

    deallocateData(data);
    free(result);

    return 0;
}
