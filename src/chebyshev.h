//
// Created by morfeush22 on 15.12.18.
//

#include "base.h"

double maxElement(struct Data data) {
    int i, j;
    double maxElement = 0;

    for (i = 0; i < data.dimension; ++i) {
        for (j = 0 ; j < data.dimension; ++j) {
            if (data.matrix[i][j] > maxElement)
                maxElement = data.matrix[i][j];
        }
    }

    return maxElement;
}

double * assignVector(double * to, const double * from, int size) {
    for (int i = 0; i < size; ++i)
        to[i] = from[i];

    return to;
}

double * addVectors(const double *from, const double *to, double * sink, int size) {
    for (int i = 0; i < size; ++i)
        sink[i] = to[i] + from[i];

    return sink;
}

double * multiplyVectorByScalar(const double * vector, const double scalar, double * sink, int size) {
    for (int i = 0; i < size; ++i)
        sink[i] = vector[i] * scalar;

    return sink;
}

double * subtractVectors(const double *from, const double *vector, double * sink, int size) {
    for (int i = 0; i < size; ++i)
        sink[i] = from[i] - vector[i];

    return sink;
}

void zeroVector(double * vector, int size) {
    for (int i = 0; i < size; ++i)
        vector[i] = 0;
}

void printVector(double * vector, int size) {
    for (int i = 0; i < size; ++i) {
        printf("[%.6f]\n", vector[i]);
    }

    printf("\n");
}

double * multiplyMatrixByVector(double ** matrix, const double * vector, double *sink, int dimension) {
    double sum;

    for (int i = 0; i < dimension; ++i) {
        sum = 0;

        for (int j = 0; j < dimension; ++j)
            sum += matrix[i][j] * vector[j];

        sink[i] = sum;
    }

    return sink;
}

double * solveLinear(struct Data data, int sParameter, int iterations) {
    int dimension = data.dimension;
    double ** matrix = data.matrix;
    double * bVector = data.bVector;

    double alfa = 100;
    double beta = 2.0 * maxElement(data);

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

    while (1) {
        printVector(xIVector, dimension);

        k = 0;
        xIVector = assignVector(xIVector, xZeroVector, dimension);
        double omegaPrev = 0;
        double omegaI = omegaZero;

        while (1) {
            t1Vector = multiplyMatrixByVector(matrix, xIVector, t1Vector, dimension);
            t1Vector = subtractVectors(t1Vector, bVector, t1Vector, dimension);
            t1Vector = multiplyVectorByScalar(t1Vector, c * (1 + omegaI * omegaPrev), t1Vector, dimension);

            t2Vector = subtractVectors(xIVector, xPrevVector, t2Vector, dimension);
            t2Vector = multiplyVectorByScalar(t2Vector, omegaI * omegaPrev, t2Vector, dimension);

            t2Vector = addVectors(xIVector, t2Vector, t2Vector, dimension);
            xPrevVector = assignVector(xPrevVector, xIVector, dimension);
            xIVector = subtractVectors(t2Vector, t1Vector, xIVector, dimension);

            omegaPrev = omegaI;
            omegaI = 1.0 / (L - omegaI);

            ++i, ++k;
            if (k >= sParameter) {
                break;
            }
        }

        xZeroVector = assignVector(xZeroVector, xPrevVector, dimension);

        if (fullIterations > iterations) {
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
