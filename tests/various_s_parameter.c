//
// Created by morfeush22 on 22.12.18.
//

#include "mpi.h"
#include "test_base.h"
#include "../src/base.h"
#include "../src/chebyshev.h"

void testDifferentSParameter() {
    int rank, size;

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    FILE *fp;
    int dimension;

    if (rank == 0) {
        fp = fopen("sources/data", "r");
        if (fp == NULL) {
            exit(EXIT_FAILURE);
        }

        dimension = numberOfLines(fp);
    }

    MPI_Bcast(&dimension, 1, MPI_INT, 0, MPI_COMM_WORLD);

    struct Data data = allocateData(dimension);
    data.dimension = dimension;

    if (rank == 0) {
        parseData(fp, &data);
        fclose(fp);
    }

    int iterations = 0;

    double * expectedResult = NULL;

    if (rank == 0) {
        expectedResult = loadVectorFromFile("sources/data.expected_solution");
    }

    for (int sParameter = 1; sParameter < 100; ++sParameter) {
        double *result = solveLinear(data, PRECISION, sParameter, &iterations, rank, size);

        if (rank == 0) {
            for (int i = 0; i < dimension; ++i) {
                ASSERT_EQUAL(expectedResult[i], result[i], PRECISION);
            }
        }

        free(result);
    }

    deallocateData(&data);
    free(expectedResult);

    MPI_Finalize();
}

int main(int argc, char ** argv) {
    testDifferentSParameter();
    printf("\n");
}