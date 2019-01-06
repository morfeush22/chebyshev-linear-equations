//
// Created by morfeush22 on 16.12.18.
//

#include "chebyshev.h"
#include "vector_operations.h"
#include "math.h"
#include "mpi.h"
#include "stdbool.h"
#include "stdlib.h"

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

double * solveLinear(const struct Data data, double precision, int sParameter, int * iterations, int rank, int size) {
    double * matrix = *data.matrix;
    double * bVector = data.bVector;
    int dimension = data.dimension;

    int chunkSize = getChunkSize(dimension, size);
    double * localMatrix = malloc(chunkSize * dimension * sizeof(double));
    double * localBVector = malloc(chunkSize * sizeof(double));

    int * counts = getCounts(dimension, size);
    int * displacements = getDisplacements(counts, size);

    int * matrixCounts = malloc(size * sizeof(int));
    int * matrixDisplacements = malloc(size * sizeof(int));

    for (int i = 0; i < size; ++i) {
        matrixCounts[i] = counts[i] * dimension;
        matrixDisplacements[i] = displacements[i] * dimension;
    }

    double * xIVector, * xZeroVector, * xPrevVector, * t1Vector, * t2Vector, * t3Vector;

    xIVector = malloc(dimension * sizeof(double));
    xZeroVector = malloc(dimension * sizeof(double));
    xPrevVector = malloc(dimension * sizeof(double));

    t1Vector = malloc(dimension * sizeof(double));
    t2Vector = malloc(dimension * sizeof(double));
    t3Vector = malloc(dimension * sizeof(double));

    MPI_Scatterv(matrix, matrixCounts, matrixDisplacements, MPI_DOUBLE, localMatrix, matrixCounts[rank], MPI_DOUBLE, 0,
                 MPI_COMM_WORLD);
    MPI_Scatterv(bVector, counts, displacements, MPI_DOUBLE, localBVector, counts[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double globalMaxMatrixElem;

    double maxMatrixElem = findMaxElementInMatrix(localMatrix, matrixCounts[rank]);
    MPI_Allreduce(&maxMatrixElem, &globalMaxMatrixElem, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    double alfa, beta, omegaZero, c, L;

    alfa = 100;
    beta = 2.0 * globalMaxMatrixElem;

    zeroVector(xIVector, dimension);
    zeroVector(xZeroVector, dimension);

    omegaZero = (beta - alfa) / (beta + alfa);
    c = 2.0 / (beta + alfa);
    L = 2.0 * (beta + alfa) / (beta - alfa);

    zeroVector(xPrevVector, dimension);

    int fullIterations = 0;

    while (true) {
        int k = 0;

        assignVector(xIVector, xZeroVector, dimension);

        double omegaPrev, omegaI;

        omegaPrev = 0;
        omegaI = omegaZero;

        for (; k < sParameter; ++k) {
            multiplyMatrixByVector(localMatrix, xIVector, t1Vector, matrixCounts[rank], dimension);
            subtractVectors(t1Vector, localBVector, t1Vector, counts[rank]);
            multiplyVectorByScalar(t1Vector, c * (1 + omegaI * omegaPrev), t1Vector, counts[rank]);

            MPI_Scatterv(xIVector, counts, displacements, MPI_DOUBLE, t3Vector, counts[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

            subtractVectors(t3Vector, xPrevVector, t2Vector, counts[rank]);
            multiplyVectorByScalar(t2Vector, omegaI * omegaPrev, t2Vector, counts[rank]);

            addVectors(t3Vector, t2Vector, t2Vector, counts[rank]);

            assignVector(xPrevVector, t3Vector, counts[rank]);
            subtractVectors(t2Vector, t1Vector, t3Vector, counts[rank]);

            omegaPrev = omegaI;
            omegaI = 1.0 / (L - omegaI);

            MPI_Allgatherv(t3Vector, counts[rank], MPI_DOUBLE, xIVector, counts, displacements, MPI_DOUBLE, MPI_COMM_WORLD);
        }

        assignVector(xZeroVector, xIVector, dimension);

        // calculate error
        multiplyMatrixByVector(localMatrix, xZeroVector, t1Vector, matrixCounts[rank], dimension);
        subtractVectors(localBVector, t1Vector, t1Vector, counts[rank]);

        double globalCurrPrecision;

        double currPrecision = findAbsMaxElementInVector(t1Vector, counts[rank]);
        MPI_Allreduce(&currPrecision, &globalCurrPrecision, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        if (globalCurrPrecision <= precision) {
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
