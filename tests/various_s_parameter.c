//
// Created by morfeush22 on 22.12.18.
//

#include "test_base.h"
#include "../src/base.h"
#include "../src/chebyshev.h"

void testDifferentSParameter() {
    struct Data data = loadDataFromFile("sources/data");
    int dimension = data.dimension;

    int iterations = 0;

    double *expectedResult = loadVectorFromFile("sources/data.expected_solution");

    for (int sParameter = 1; sParameter < 100; ++sParameter) {
        double *result = solveLinear(data, PRECISION, sParameter, &iterations);

        for (int i = 0; i < dimension; ++i) {
            ASSERT_EQUAL(expectedResult[i], result[i], PRECISION);
        }

        free(result);
    }

    deallocateData(data);
    free(expectedResult);
}

int main(int argc, char ** argv) {
    testDifferentSParameter();
    printf("\n");
}