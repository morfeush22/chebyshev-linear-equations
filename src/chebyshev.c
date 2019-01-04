//
// Created by morfeush22 on 16.12.18.
//

#include "chebyshev.h"
#include "vector_operations.h"
#include "math.h"
#include "mpi.h"
#include "stdbool.h"
#include "stdlib.h"

double * solveLinear(const struct Data data, double precision, int sParameter, int * iterations, int rank, int size) {
    const double * const * matrix = (const double * const *)data.matrix;
    const double * bVector = data.bVector;
    int dimension = data.dimension;

    double alfa = 100;
    double beta = 2.0 * findMaxElementInMatrix(matrix, dimension);

    double * xIVector, * xZeroVector, * xPrevVector, * t1Vector, * t2Vector;
    xIVector = xZeroVector = xPrevVector = t1Vector = t2Vector = NULL;

    double omegaZero, c, L;

    if (rank == 0) {
        xIVector = malloc(dimension * sizeof(double));
        xZeroVector = malloc(dimension * sizeof(double));
        xPrevVector = malloc(dimension * sizeof(double));
        zeroVector(xIVector, dimension);
        zeroVector(xZeroVector, dimension);

        t1Vector = malloc(dimension * sizeof(double));
        t2Vector = malloc(dimension * sizeof(double));

        omegaZero = (beta - alfa) / (beta + alfa);
        c = 2.0 / (beta + alfa);
        L = 2.0 * (beta + alfa) / (beta - alfa);

        zeroVector(xPrevVector, dimension);
    }

    int i = 0;
    int fullIterations = 0;

    while (true) {
        int k = 0;

        assignVector(xIVector, xZeroVector, dimension);

        double omegaPrev, omegaI;

        if (rank == 0) {
            omegaPrev = 0;
            omegaI = omegaZero;
        }

        for (; k < sParameter; ++i, ++k) {
            multiplyMatrixByVector(matrix, xIVector, t1Vector, dimension);
            subtractVectors(t1Vector, bVector, t1Vector, dimension);
            multiplyVectorByScalar(t1Vector, c * (1 + omegaI * omegaPrev), t1Vector, dimension);

            subtractVectors(xIVector, xPrevVector, t2Vector, dimension);
            multiplyVectorByScalar(t2Vector, omegaI * omegaPrev, t2Vector, dimension);

            addVectors(xIVector, t2Vector, t2Vector, dimension);

            assignVector(xPrevVector, xIVector, dimension);
            subtractVectors(t2Vector, t1Vector, xIVector, dimension);

            if (rank == 0) {
                omegaPrev = omegaI;
                omegaI = 1.0 / (L - omegaI);
            }
        }

        assignVector(xZeroVector, xIVector, dimension);

        // calculate error
        multiplyMatrixByVector(matrix, xZeroVector, t1Vector, dimension);
        subtractVectors(bVector, t1Vector, t1Vector, dimension);

        double currPrecision = fabs(findAbsMaxElementInVector(t1Vector, dimension));

        MPI_Bcast(&currPrecision, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        if (currPrecision <= precision) {
            *iterations = fullIterations;
            break;
        }

        ++fullIterations;
    }

    free(xZeroVector);
    free(xPrevVector);
    free(t1Vector);
    free(t2Vector);

    if (rank == 0) {
        return xIVector;
    }
    else {
        return NULL;
    }
}
