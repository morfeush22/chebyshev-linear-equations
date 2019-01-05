//
// Created by morfeush22 on 16.12.18.
//

#include "base.h"
#include "stdlib.h"

void deallocateData(struct Data * data) {
    double * memBlock = data->matrix[0];

    free(memBlock);
    free(data->matrix);
    free(data->bVector);
}

void printData(struct Data * data) {
    int i, j;

    for (i = 0; i < data->dimension; ++i) {
        for (j = 0 ; j < data->dimension; ++j) {
            printf("%.6f ", data->matrix[i][j]);
        }

        printf("%.6f\n", data->bVector[i]);
    }
}

void parseData(FILE * fp, struct Data * data) {
    readData(fp, data);
}

struct Data allocateData(int dimension) {
    double * memBlock = malloc(dimension * dimension * sizeof(double));
    double ** matrix = malloc(dimension * sizeof(double *));

    for (int i = 0; i < dimension; ++i) {
        matrix[i] = memBlock + i * dimension;
    }

    double * bVector = malloc(dimension * sizeof(double));

    struct Data data;
    data.matrix = matrix;
    data.bVector = bVector;
    data.dimension = dimension;

    return data;
}

void readData(FILE * fp, struct Data * data) {
    int i, j;
    double var;

    for (i = 0; i < data->dimension; ++i) {
        for (j = 0 ; j < data->dimension; ++j) {
            fscanf(fp, "%lf", &var);
            data->matrix[i][j] = var;
        }

        fscanf(fp,"%lf", &var);
        data->bVector[i] = var;
    }
}

int numberOfLines(FILE * fp) {
    int ch;
    int lines = 0;

    while(!feof(fp)) {
        ch = fgetc(fp);

        if (ch == '\n') {
            lines++;
        }
    }

    fseek(fp, 0, SEEK_SET);

    return lines;
}
