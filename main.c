#include "src/chebyshev.h"

int main() {
    int sParameter = 100;
    double iterations = 20;

    double * result = solveLinear(loadDataFromFile("tests/sources/simple_eq"), sParameter, iterations);

    return 0;
}