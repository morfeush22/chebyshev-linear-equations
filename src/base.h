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

struct Data loadDataFromFile(const char * path);
void deallocateData(struct Data data);

double maxMatrixElement(struct Data data);

void printData(struct Data data);


static struct Data parseData(FILE * fp, int dimension);
static struct Data allocateData(int dimension);
static void readData(FILE * fp, struct Data data);
static int numberOfLines(FILE * fp);

#endif //CHEBYSHEV_BASE_H
