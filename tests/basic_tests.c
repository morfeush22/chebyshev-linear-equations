//
// Created by morfeush22 on 15.12.18.
//

#include "test_base.h"
#include "../src/base.h"
#include "../src/chebyshev.h"

void basicTest() {
    struct Data data = loadDataFromFile("sources/data");
    int dimension = data.dimension;

    int sParameter = 10;
    int iterations = 0;

    double * result = solveLinear(data, PRECISION, sParameter, &iterations);
    double * expectedResult = loadVectorFromFile("sources/data.expected_solution");

    for (int i = 0; i < dimension; ++i) {
        ASSERT_EQUAL(expectedResult[i], result[i], PRECISION);
    }

    deallocateData(data);
    free(result);
    free(expectedResult);
}

void testDifferentSParameter() {
    struct Data data = loadDataFromFile("sources/data");
    int dimension = data.dimension;

    int iterations = 0;

    double *expectedResult = loadVectorFromFile("sources/data.expected_solution");

    for (int sParameter = 2; sParameter < 1000; ++sParameter) {
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
    basicTest();
    testDifferentSParameter();
    printf("\n");
}
