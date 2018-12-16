//
// Created by morfeush22 on 16.12.18.
//

#ifndef CHEBYSHEV_VECTOR_OPERATIONS_H
#define CHEBYSHEV_VECTOR_OPERATIONS_H

void addVectors(const double *from, const double *to, double * sink, int size);
void assignVector(double * to, const double * from, int size);
void multiplyMatrixByVector(double ** matrix, const double * vector, double *sink, int dimension);
void multiplyVectorByScalar(const double * vector, double scalar, double * sink, int size);
void printVector(double * vector, int size);
void subtractVectors(const double *from, const double *vector, double * sink, int size);
void zeroVector(double * vector, int size);

#endif //CHEBYSHEV_VECTOR_OPERATIONS_H
