//
// Created by morfeush22 on 16.12.18.
//

#include "vector_operations.h"
#include "stdio.h"

void addVectors(const double *from, const double *to, double * sink, int size) {
    for (int i = 0; i < size; ++i)
        sink[i] = to[i] + from[i];
}

void assignVector(double * to, const double * from, int size) {
    for (int i = 0; i < size; ++i)
        to[i] = from[i];
}

void multiplyMatrixByVector(double ** matrix, const double * vector, double *sink, int dimension) {
    double sum;

    for (int i = 0; i < dimension; ++i) {
        sum = 0;

        for (int j = 0; j < dimension; ++j)
            sum += matrix[i][j] * vector[j];

        sink[i] = sum;
    }
}

void multiplyVectorByScalar(const double * vector, double scalar, double * sink, int size) {
    for (int i = 0; i < size; ++i)
        sink[i] = vector[i] * scalar;
}

void printVector(double * vector, int size) {
    for (int i = 0; i < size; ++i) {
        printf("[%.6f]\n", vector[i]);
    }

    printf("\n");
}

void subtractVectors(const double *from, const double *vector, double * sink, int size) {
    for (int i = 0; i < size; ++i)
        sink[i] = from[i] - vector[i];
}

void zeroVector(double * vector, int size) {
    for (int i = 0; i < size; ++i)
        vector[i] = 0;
}
