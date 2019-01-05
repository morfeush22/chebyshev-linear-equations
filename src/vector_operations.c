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
    if (rank == 0) {
        for (int i = 0; i < size; ++i) {
            sink[i] = vector1[i] + vector2[i];
        }
    }
}

void assignVector(double * to, const double * from, int size) {
    for (int i = 0; i < size; ++i) {
        to[i] = from[i];
    }
}

double findMaxElementInMatrix(const double * const * matrix, int dimension, int rank, int procNum) {
    if (rank == 0) {
        double maxElement = 0;

        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                maxElement = fmax(maxElement, matrix[i][j]);
            }
        }

        return maxElement;
    }
}

double findAbsMaxElementInVector(const double * vector, int size, int rank, int procNum) {
    if (rank == 0) {
        double maxElement = 0;

        for (int i = 0; i < size; ++i) {
            maxElement = fmax(maxElement, fabs(vector[i]));
        }

        return maxElement;
    }
}

void multiplyMatrixByVector(const double * const * matrix, const double * vector, double * sink, int dimension,
        int rank, int procNum) {
    if (rank == 0) {
        for (int i = 0; i < dimension; ++i) {
            double sum = 0;

            for (int j = 0; j < dimension; ++j) {
                sum += matrix[i][j] * vector[j];
            }

            sink[i] = sum;
        }
    }
}

void multiplyVectorByScalar(const double * vector, double scalar, double * sink, int size, int rank, int procNum) {
    if (rank == 0) {
        for (int i = 0; i < size; ++i) {
            sink[i] = vector[i] * scalar;
        }
    }
}

void printVector(double * vector, int size) {
    for (int i = 0; i < size; ++i) {
        printf("%.6f\n", vector[i]);
    }
}

void subtractVectors(const double * from, const double * vector, double * sink, int size, int rank, int procNum) {
    MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);

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
