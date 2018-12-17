//
// Created by morfeush22 on 16.12.18.
//

#include "vector_operations.h"
#include "stdio.h"

void addVectors(const double * vector1, const double * vector2, double * sink, int size) {
    for (int i = 0; i < size; ++i) {
        sink[i] = vector1[i] + vector2[i];
    }
}

void assignVector(double * to, const double * from, int size) {
    for (int i = 0; i < size; ++i) {
        to[i] = from[i];
    }
}

double findMaxElementInMatrix(const double * const * matrix, int dimension) {
    int i, j;
    double maxElement = 0;

    for (i = 0; i < dimension; ++i) {
        for (j = 0 ; j < dimension; ++j) {
            if (matrix[i][j] > maxElement) {
                maxElement = matrix[i][j];
            }
        }
    }

    return maxElement;
}

double findMaxElementInVector(const double * vector, int size) {
    double maxElement = 0;

    for (int i = 0; i < size; ++i) {
        if (vector[i] > maxElement) {
            maxElement = vector[i];
        }
    }

    return maxElement;
}

void multiplyMatrixByVector(const double * const* matrix, const double * vector, double *sink, int dimension) {
    double sum;

    for (int i = 0; i < dimension; ++i) {
        sum = 0;

        for (int j = 0; j < dimension; ++j)
            sum += matrix[i][j] * vector[j];

        sink[i] = sum;
    }
}

void multiplyVectorByScalar(const double * vector, double scalar, double * sink, int size) {
    for (int i = 0; i < size; ++i) {
        sink[i] = vector[i] * scalar;
    }
}

void printVector(double * vector, int size) {
    for (int i = 0; i < size; ++i) {
        printf("[%.6f]\n", vector[i]);
    }

    printf("\n");
}

void subtractVectors(const double * from, const double * vector, double * sink, int size) {
    for (int i = 0; i < size; ++i) {
        sink[i] = from[i] - vector[i];
    }
}

void zeroVector(double * vector, int size) {
    for (int i = 0; i < size; ++i) {
        vector[i] = 0;
    }
}
