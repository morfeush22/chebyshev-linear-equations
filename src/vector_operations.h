//
// Created by morfeush22 on 16.12.18.
//

#ifndef CHEBYSHEV_VECTOR_OPERATIONS_H
#define CHEBYSHEV_VECTOR_OPERATIONS_H

void addVectors(const double * vector1, const double * vector2, double * sink, int size);
void assignVector(double * to, const double * from, int size);
double findMaxElementInMatrix(const double * const * matrix, int dimension);
double findAbsMaxElementInVector(const double * vector, int size);
void multiplyMatrixByVector(const double * const * matrix, const double * vector, double * sink, int dimension);
void multiplyVectorByScalar(const double * vector, double scalar, double * sink, int size);
void printVector(double * vector, int size);
void subtractVectors(const double * from, const double * vector, double * sink, int size);
void zeroVector(double * vector, int size);

#endif //CHEBYSHEV_VECTOR_OPERATIONS_H
