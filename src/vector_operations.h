//
// Created by morfeush22 on 16.12.18.
//

#ifndef CHEBYSHEV_VECTOR_OPERATIONS_H
#define CHEBYSHEV_VECTOR_OPERATIONS_H

void addVectors(const double * vector1, const double * vector2, double * sink, int size, int rank, int procNum);
void assignVector(double * to, const double * from, int size);
double findMaxElementInMatrix(const double * const * matrix, int dimension, int rank, int procNum);
double findAbsMaxElementInVector(const double * vector, int size, int rank, int procNum);
void multiplyMatrixByVector(const double * const * matrix, double * vector, double * sink, int dimension,
        int rank, int procNum);
void multiplyVectorByScalar(const double * vector, double scalar, double * sink, int size, int rank, int procNum);
void printVector(double * vector, int size);
void subtractVectors(const double * from, const double * vector, double * sink, int size, int rank, int procNum);
void zeroVector(double * vector, int size);

#endif //CHEBYSHEV_VECTOR_OPERATIONS_H
