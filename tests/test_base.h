//
// Created by morfeush22 on 15.12.18.
//

#ifndef CHEBYSHEV_TEST_BASE_H
#define CHEBYSHEV_TEST_BASE_H

#include "stdlib.h"
#include "stdio.h"

#define ASSERT_EQUAL(x, y)                              \
{                                                       \
    if ( (x) != (y) )                                   \
    {                                                   \
        printf("%.6f IS DIFFERENT FROM %.6f\n", y, x);  \
        exit(EXIT_FAILURE);                             \
    }                                                   \
}

#endif //CHEBYSHEV_TEST_BASE_H
