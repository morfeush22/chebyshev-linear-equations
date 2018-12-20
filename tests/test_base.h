//
// Created by morfeush22 on 15.12.18.
//

#ifndef CHEBYSHEV_TEST_BASE_H
#define CHEBYSHEV_TEST_BASE_H

#include "../src/base.h"
#include "math.h"
#include "stdlib.h"
#include "stdio.h"

#define PRECISION 0.000001

#define ASSERT_EQUAL(x, y, precision)                                           \
{                                                                               \
    if ( fabs((x) - (y)) > (precision))                                         \
    {                                                                           \
        printf("%.6f IS DIFFERENT FROM %.6f more than %.6f\n", y, x, precision);\
        exit(EXIT_FAILURE);                                                     \
    }                                                                           \
    printf(".");                                                                \
}

static void readVector(FILE * fp, double * sink, int size) {
    double var;

    for (int i = 0; i < size; ++i) {
        fscanf(fp, "%lf", &var);
        sink[i] = var;
    }
}

static double * parseVector(FILE * fp, int size) {
    double * vector = malloc(size * sizeof(double));

    readVector(fp, vector, size);

    return vector;
}

double * loadVectorFromFile(const char * path) {
    FILE * fp;

    fp = fopen(path, "r");
    if (fp == NULL) {
        exit(EXIT_FAILURE);
    }

    int size = numberOfLines(fp);

    double * vector = parseVector(fp, size);

    fclose(fp);

    return vector;
}

#endif //CHEBYSHEV_TEST_BASE_H
