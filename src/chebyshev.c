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

    double * xIVector, * xZeroVector, * xPrevVector, * t1Vector, * t2Vector;

    xIVector = malloc(dimension * sizeof(double));
    xZeroVector = malloc(dimension * sizeof(double));
    xPrevVector = malloc(dimension * sizeof(double));

    t1Vector = malloc(dimension * sizeof(double));
    t2Vector = malloc(dimension * sizeof(double));

    double maxMatrixElem = findMaxElementInMatrix(matrix, dimension, rank, size);
    double alfa, beta, omegaZero, c, L;

    if (rank == 0) {
        alfa = 100;
        beta = 2.0 * maxMatrixElem;

        zeroVector(xIVector, dimension);
        zeroVector(xZeroVector, dimension);

        omegaZero = (beta - alfa) / (beta + alfa);
        c = 2.0 / (beta + alfa);
        L = 2.0 * (beta + alfa) / (beta - alfa);

        zeroVector(xPrevVector, dimension);
    }

    int fullIterations = 0;

    while (true) {
        int k = 0;

        if (rank == 0) {
            assignVector(xIVector, xZeroVector, dimension);
        }

        double omegaPrev, omegaI;

        if (rank == 0) {
            omegaPrev = 0;
            omegaI = omegaZero;
        }

        for (; k < sParameter; ++k) {
            multiplyMatrixByVector(matrix, xIVector, t1Vector, dimension, rank, size);
            subtractVectors(t1Vector, bVector, t1Vector, dimension, rank, size);
            multiplyVectorByScalar(t1Vector, c * (1 + omegaI * omegaPrev), t1Vector, dimension, rank, size);

            subtractVectors(xIVector, xPrevVector, t2Vector, dimension, rank, size);
            multiplyVectorByScalar(t2Vector, omegaI * omegaPrev, t2Vector, dimension, rank, size);

            addVectors(xIVector, t2Vector, t2Vector, dimension, rank, size);

            assignVector(xPrevVector, xIVector, dimension);
            subtractVectors(t2Vector, t1Vector, xIVector, dimension, rank, size);

            if (rank == 0) {
                omegaPrev = omegaI;
                omegaI = 1.0 / (L - omegaI);
            }
        }

        if (rank == 0) {
            assignVector(xZeroVector, xIVector, dimension);
        }

        // calculate error
        multiplyMatrixByVector(matrix, xZeroVector, t1Vector, dimension, rank, size);
        subtractVectors(bVector, t1Vector, t1Vector, dimension, rank, size);

        double currPrecision = findAbsMaxElementInVector(t1Vector, dimension, rank, size);

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

    return xIVector;
}
