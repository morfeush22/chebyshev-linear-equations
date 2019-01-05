//
// Created by morfeush22 on 15.12.18.
//

#ifndef CHEBYSHEV_BASE_H
#define CHEBYSHEV_BASE_H

#include "stdio.h"

struct Data {
    double ** matrix;
    double * bVector;
    int dimension;
};

struct Data allocateData(int dimension);
void deallocateData(struct Data * data);
void parseData(FILE * fp, struct Data * data);
void printData(struct Data * data);

int numberOfLines(FILE * fp);

static void readData(FILE * fp, struct Data * data);

#endif //CHEBYSHEV_BASE_H
