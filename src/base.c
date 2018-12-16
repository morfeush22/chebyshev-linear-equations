//
// Created by morfeush22 on 16.12.18.
//

#include "base.h"
#include "stdlib.h"

struct Data loadDataFromFile(const char * path) {
    FILE * fp;

    fp = fopen(path, "r");
    if (fp == NULL)
        exit(EXIT_FAILURE);

    int dimension = numberOfLines(fp);

    struct Data data = parseData(fp, dimension);

    fclose(fp);

    return data;
}

void deallocateData(struct Data data) {
    for (int i = 0; i < data.dimension; ++i)
        free(data.matrix[i]);

    free(data.matrix);
    free(data.bVector);
}

double maxMatrixElement(struct Data data) {
    int i, j;
    double maxElement = 0;

    for (i = 0; i < data.dimension; ++i) {
        for (j = 0 ; j < data.dimension; ++j) {
            if (data.matrix[i][j] > maxElement) {
                maxElement = data.matrix[i][j];
            }
        }
    }

    return maxElement;
}

void printData(struct Data data) {
    int i, j;

    for (i = 0; i < data.dimension; ++i) {
        printf("[ ");
        for (j = 0 ; j < data.dimension; ++j) {
            printf("%.6f ", data.matrix[i][j]);
        }

        printf("][ ");
        printf("%.6f ", data.bVector[i]);
        printf("]\n");
    }
}

struct Data parseData(FILE * fp, int dimension) {
    struct Data data = allocateData(dimension);

    readData(fp, data);

    return data;
}

struct Data allocateData(int dimension) {
    double ** matrix = malloc(dimension * sizeof(double *));

    for (int i = 0; i < dimension; ++i) {
        matrix[i] = malloc(dimension * sizeof(double));
    }

    double * bVector = malloc(dimension * sizeof(double));

    struct Data data;
    data.matrix = matrix;
    data.bVector = bVector;
    data.dimension = dimension;

    return data;
}

void readData(FILE * fp, struct Data data) {
    int i, j;
    double var;

    for (i = 0; i < data.dimension; ++i) {
        for (j = 0 ; j < data.dimension; ++j) {
            fscanf(fp,"%lf", &var);
            data.matrix[i][j] = var;
        }

        fscanf(fp,"%lf",&var);
        data.bVector[i] = var;
    }
}

int numberOfLines(FILE * fp) {
    int ch;
    int lines = 0;

    while(!feof(fp)) {
        ch = fgetc(fp);

        if (ch == '\n')
            lines++;
    }

    fseek(fp, 0, SEEK_SET);

    return lines;
}
