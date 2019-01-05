//
// Created by morfeush22 on 16.12.18.
//

#include "vector_operations.h"
#include "math.h"
#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"

int getChunkSize(int size, int procNum) {
    return size % procNum ? size / procNum + 1 : size / procNum;
}

int * getCounts(int size, int procNum) {
    int * counts = malloc(procNum * sizeof(int));

    for (int i = 0; i < procNum; ++i) {
        counts[i] = getChunkSize(size, procNum - i);
        size -= counts[i];
    }

    return counts;
}

int * getDisplacements(const int * counts, int procNum) {
    int * displacements = malloc(procNum * sizeof(int));

    int cum = 0;

    for (int i = 0; i < procNum; ++i) {
        displacements[i] = cum;
        cum += counts[i];
    }

    return displacements;
}

void addVectors(const double * vector1, const double * vector2, double * sink, int size, int rank, int procNum) {
    size_t mallocSize = getChunkSize(size, procNum) * sizeof(double);
    double * localVec1 = malloc(mallocSize);
    double * localVec2 = malloc(mallocSize);
    double * localSink = malloc(mallocSize);

    int * counts = getCounts(size, procNum);
    int * displacements = getDisplacements(counts, procNum);

    MPI_Scatterv(vector1, counts, displacements, MPI_DOUBLE, localVec1, counts[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(vector2, counts, displacements, MPI_DOUBLE, localVec2, counts[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for (int i = 0; i < counts[rank]; ++i) {
        localSink[i] = localVec1[i] + localVec2[i];
    }

    MPI_Gatherv(localSink, counts[rank], MPI_DOUBLE, sink, counts, displacements, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    free(localVec1);
    free(localVec2);
    free(localSink);

    free(counts);
    free(displacements);
}

void assignVector(double * to, const double * from, int size) {
    for (int i = 0; i < size; ++i) {
        to[i] = from[i];
    }
}

double findMaxElementInMatrix(const double * const * matrix, int dimension, int rank, int procNum) {
    size_t mallocSize = getChunkSize(dimension, procNum) * sizeof(double);
    double * localMatrix = malloc(mallocSize * dimension);

    int * counts = getCounts(dimension, procNum);
    int * displacements = getDisplacements(counts, procNum);

    int * matrixCounts = malloc(procNum * sizeof(int));
    int * matrixDisplacements = malloc(procNum * sizeof(int));

    for (int i = 0; i < procNum; ++i) {
        matrixCounts[i] = counts[i] * dimension;
        matrixDisplacements[i] = displacements[i] * dimension;
    }

    MPI_Scatterv(*matrix, matrixCounts, matrixDisplacements, MPI_DOUBLE, localMatrix, matrixCounts[rank], MPI_DOUBLE, 0,
                 MPI_COMM_WORLD);

    double maxElement = 0;

    for (int i = 0; i < counts[rank]; ++i) {
        for (int j = 0; j < dimension; ++j) {
            maxElement = fmax(maxElement, *(localMatrix + i * dimension + j));
        }
    }

    double globalMaxElement;

    MPI_Reduce(&maxElement, &globalMaxElement, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    free(localMatrix);

    free(counts);
    free(displacements);

    free(matrixCounts);
    free(matrixDisplacements);

    return globalMaxElement;
}

double findAbsMaxElementInVector(const double * vector, int size, int rank, int procNum) {
    size_t mallocSize = getChunkSize(size, procNum) * sizeof(double);
    double * localVector = malloc(mallocSize);

    int * counts = getCounts(size, procNum);
    int * displacements = getDisplacements(counts, procNum);

    MPI_Scatterv(vector, counts, displacements, MPI_DOUBLE, localVector, counts[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double maxElement = 0;

    for (int i = 0; i < counts[rank]; ++i) {
        maxElement = fmax(maxElement, fabs(localVector[i]));
    }

    double globalMaxElement;

    MPI_Reduce(&maxElement, &globalMaxElement, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    free(localVector);

    free(counts);
    free(displacements);

    return globalMaxElement;
}

void multiplyMatrixByVector(const double * const * matrix, double * vector, double * sink, int dimension, int rank,
        int procNum) {
    size_t mallocSize = getChunkSize(dimension, procNum) * sizeof(double);
    double * localMatrix = malloc(mallocSize * dimension);
    double * localSink = malloc(mallocSize);

    int * counts = getCounts(dimension, procNum);
    int * displacements = getDisplacements(counts, procNum);

    int * matrixCounts = malloc(procNum * sizeof(int));
    int * matrixDisplacements = malloc(procNum * sizeof(int));

    for (int i = 0; i < procNum; ++i) {
        matrixCounts[i] = counts[i] * dimension;
        matrixDisplacements[i] = displacements[i] * dimension;
    }

    MPI_Bcast(vector, dimension, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(*matrix, matrixCounts, matrixDisplacements, MPI_DOUBLE, localMatrix, matrixCounts[rank], MPI_DOUBLE, 0,
            MPI_COMM_WORLD);

    for (int i = 0; i < counts[rank]; ++i) {
        double sum = 0;

        for (int j = 0; j < dimension; ++j) {
            sum += *(localMatrix + i * dimension + j) * vector[j];
        }

        localSink[i] = sum;
    }

    MPI_Gatherv(localSink, counts[rank], MPI_DOUBLE, sink, counts, displacements, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    free(localMatrix);
    free(localSink);

    free(counts);
    free(displacements);

    free(matrixCounts);
    free(matrixDisplacements);
}

void multiplyVectorByScalar(const double * vector, double scalar, double * sink, int size, int rank, int procNum) {
    size_t mallocSize = getChunkSize(size, procNum) * sizeof(double);
    double * localVector = malloc(mallocSize);
    double * localSink = malloc(mallocSize);

    int * counts = getCounts(size, procNum);
    int * displacements = getDisplacements(counts, procNum);

    MPI_Bcast(&scalar, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(vector, counts, displacements, MPI_DOUBLE, localVector, counts[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for (int i = 0; i < counts[rank]; ++i) {
        localSink[i] = localVector[i] * scalar;
    }

    MPI_Gatherv(localSink, counts[rank], MPI_DOUBLE, sink, counts, displacements, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    free(localVector);
    free(localSink);

    free(counts);
    free(displacements);
}

void printVector(double * vector, int size) {
    for (int i = 0; i < size; ++i) {
        printf("%.6f\n", vector[i]);
    }
}

void subtractVectors(const double * from, const double * vector, double * sink, int size, int rank, int procNum) {
    size_t mallocSize = getChunkSize(size, procNum) * sizeof(double);
    double * localFrom = malloc(mallocSize);
    double * localVector = malloc(mallocSize);
    double * localSink = malloc(mallocSize);

    int * counts = getCounts(size, procNum);
    int * displacements = getDisplacements(counts, procNum);

    MPI_Scatterv(from, counts, displacements, MPI_DOUBLE, localFrom, counts[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(vector, counts, displacements, MPI_DOUBLE, localVector, counts[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for (int i = 0; i < counts[rank]; ++i) {
        localSink[i] = localFrom[i] - localVector[i];
    }

    MPI_Gatherv(localSink, counts[rank], MPI_DOUBLE, sink, counts, displacements, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    free(localFrom);
    free(localVector);
    free(localSink);

    free(counts);
    free(displacements);
}

void zeroVector(double * vector, int size) {
    memset(vector, 0, size * sizeof(double));
}
