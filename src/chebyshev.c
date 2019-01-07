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

double * solveLinear(const struct Data * data, double precision, int sParameter, int * iterations, int rank, int size) {
    double * matrix = *data->matrix;
    double * bVector = data->bVector;
    int dimension = data->dimension;

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

    MPI_Scatterv(matrix, matrixCounts, matrixDisplacements, MPI_DOUBLE, localMatrix, matrixCounts[rank], MPI_DOUBLE, 0,
                 MPI_COMM_WORLD);
    MPI_Scatterv(bVector, counts, displacements, MPI_DOUBLE, localBVector, counts[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double * xIVector = malloc(dimension * sizeof(double));
    double * xZeroVector = malloc(dimension * sizeof(double));

    zeroVector(xIVector, dimension);
    zeroVector(xZeroVector, dimension);

    double * localXIVector = malloc(chunkSize * sizeof(double));
    double * xPrevVector = malloc(chunkSize * sizeof(double));
    double * t1Vector = malloc(chunkSize * sizeof(double));
    double * t2Vector = malloc(chunkSize * sizeof(double));

    zeroVector(xPrevVector, chunkSize);

    double globalMaxMatrixElem;

    double maxMatrixElem = findMaxElementInMatrix(localMatrix, matrixCounts[rank]);
    MPI_Allreduce(&maxMatrixElem, &globalMaxMatrixElem, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    double alfa = 100;
    double beta = 2.0 * globalMaxMatrixElem;

    double omegaZero = (beta - alfa) / (beta + alfa);
    double c = 2.0 / (beta + alfa);
    double L = 2.0 * (beta + alfa) / (beta - alfa);

    int fullIterations = 0;

    while (true) {
        int k = 0;

        assignVector(xIVector, xZeroVector, dimension);

        double omegaPrev = 0;
        double omegaI = omegaZero;

        for (; k < sParameter; ++k) {
            multiplyMatrixByVector(localMatrix, xIVector, t1Vector, matrixCounts[rank], dimension);
            subtractVectors(t1Vector, localBVector, t1Vector, counts[rank]);
            multiplyVectorByScalar(t1Vector, c * (1 + omegaI * omegaPrev), t1Vector, counts[rank]);

            MPI_Scatterv(xIVector, counts, displacements, MPI_DOUBLE, localXIVector, counts[rank], MPI_DOUBLE, 0,
                    MPI_COMM_WORLD);

            subtractVectors(localXIVector, xPrevVector, t2Vector, counts[rank]);
            multiplyVectorByScalar(t2Vector, omegaI * omegaPrev, t2Vector, counts[rank]);

            addVectors(localXIVector, t2Vector, t2Vector, counts[rank]);

            assignVector(xPrevVector, localXIVector, counts[rank]);
            subtractVectors(t2Vector, t1Vector, localXIVector, counts[rank]);

            omegaPrev = omegaI;
            omegaI = 1.0 / (L - omegaI);

            MPI_Allgatherv(localXIVector, counts[rank], MPI_DOUBLE, xIVector, counts, displacements, MPI_DOUBLE,
                    MPI_COMM_WORLD);
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

    free(localMatrix);
    free(localBVector);

    free(counts);
    free(displacements);

    free(matrixCounts);
    free(matrixDisplacements);

    free(xZeroVector);

    free(localXIVector);
    free(xPrevVector);
    free(t1Vector);
    free(t2Vector);

    return xIVector;
}
