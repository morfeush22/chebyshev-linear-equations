//
// Created by morfeush22 on 16.12.18.
//

#include "vector_operations.h"
#include "math.h"
#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"

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

double findMaxElementInMatrix(const double * matrix, int size) {
    double maxElement = 0;

    for (int i = 0; i < size; ++i) {
        maxElement = fmax(maxElement, *(matrix + i));
    }

    return maxElement;
}

double findAbsMaxElementInVector(const double * vector, int size) {
    double maxElement = 0;

    for (int i = 0; i < size; ++i) {
        maxElement = fmax(maxElement, fabs(vector[i]));
    }

    return maxElement;
}

void multiplyMatrixByVector(const double * matrix, const double * vector, double * sink, int matrixSize, int vectorSize) {
    int numRows = matrixSize / vectorSize;

    for (int i = 0; i < numRows; ++i) {
        double sum = 0;

        for (int j = 0; j < vectorSize; ++j) {
            sum += *(matrix + i * vectorSize + j) * vector[j];
        }

        sink[i] = sum;
    }
}

void multiplyVectorByScalar(const double * vector, double scalar, double * sink, int size) {
    for (int i = 0; i < size; ++i) {
        sink[i] = vector[i] * scalar;
    }
}

void printVector(const double * vector, int size) {
    for (int i = 0; i < size; ++i) {
        printf("%.6f\n", vector[i]);
    }
}

void subtractVectors(const double * from, const double * vector, double * sink, int size) {
    for (int i = 0; i < size; ++i) {
        sink[i] = from[i] - vector[i];
    }
}

void zeroVector(double * vector, int size) {
    memset(vector, 0, size * sizeof(double));
}
