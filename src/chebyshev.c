//
// Created by morfeush22 on 16.12.18.
//

#include "chebyshev.h"
#include "vector_operations.h"
#include "math.h"
#include "stdbool.h"
#include "stdlib.h"

double * solveLinear(const struct Data data, double precision, int sParameter, int * iterations) {
    const double * const * matrix = (const double * const *)data.matrix;
    const double * bVector = data.bVector;
    int dimension = data.dimension;

    double alfa = 100;
    double beta = 2.0 * findMaxElementInMatrix(matrix, dimension);

    double * xIVector = malloc(dimension * sizeof(double));
    double * xZeroVector = malloc(dimension * sizeof(double));
    double * xPrevVector = malloc(dimension * sizeof(double));
    zeroVector(xIVector, dimension);
    zeroVector(xZeroVector, dimension);

    double * t1Vector = malloc(dimension * sizeof(double));
    double * t2Vector = malloc(dimension * sizeof(double));

    int i, k;

    i = 0;
    double omegaZero = (beta - alfa) / (beta + alfa);
    double c = 2.0 / (beta + alfa);
    double L = 2.0 * (beta + alfa) / (beta - alfa);
    zeroVector(xPrevVector, dimension);

    int fullIterations = 0;

    while (true) {
#ifdef DEBUG_PRINTS
        printf("X VECTOR\n");
        printVector(xIVector, dimension);
#endif

        k = 0;
        assignVector(xIVector, xZeroVector, dimension);
        double omegaPrev = 0;
        double omegaI = omegaZero;

        while (true) {
            multiplyMatrixByVector(matrix, xIVector, t1Vector, dimension);
            subtractVectors(t1Vector, bVector, t1Vector, dimension);
            multiplyVectorByScalar(t1Vector, c * (1 + omegaI * omegaPrev), t1Vector, dimension);

            subtractVectors(xIVector, xPrevVector, t2Vector, dimension);
            multiplyVectorByScalar(t2Vector, omegaI * omegaPrev, t2Vector, dimension);

            addVectors(xIVector, t2Vector, t2Vector, dimension);
            assignVector(xPrevVector, xIVector, dimension);
            subtractVectors(t2Vector, t1Vector, xIVector, dimension);

            omegaPrev = omegaI;
            omegaI = 1.0 / (L - omegaI);

            ++i, ++k;
            if (k >= sParameter) {
                break;
            }
        }

        assignVector(xZeroVector, xIVector, dimension);

        // calculate error
        multiplyMatrixByVector(matrix, xZeroVector, t1Vector, dimension);
        subtractVectors(bVector, t1Vector, t1Vector, dimension);

#ifdef DEBUG_PRINTS
        printf("DIFF VECTOR\n");
        printVector(t1Vector, dimension);
#endif

        double currPrecision = fabs(findAbsMaxElementInVector(t1Vector, dimension));

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
