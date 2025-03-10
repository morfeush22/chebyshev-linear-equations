#!/usr/bin/env python

import sys
import numpy


def format_output(value):
    return "%.6f" % value

def rvs(loc, scale, dimension):
    random_state = numpy.random

    H = numpy.eye(dimension)
    D = numpy.ones((dimension,))

    for n in range(1, dimension):
        x = random_state.normal(loc=loc, scale=scale, size=(dimension - n + 1,))
        D[n - 1] = numpy.sign(x[0])
        x[0] = x[0] - D[n - 1] * numpy.sqrt((x * x).sum())
        Hx = (numpy.eye(dimension - n + 1) - 2.0 * numpy.outer(x, x) / (x * x).sum())
        mat = numpy.eye(dimension)
        mat[n - 1:, n - 1:] = Hx
        H = numpy.dot(H, mat)

    D[-1] = (-1) ** (1 - (dimension % 2)) * D.prod()

    H = (D * H.T).T

    return H

def generate_a_matrix(loc, scale, a, b, dimension):
    mat = numpy.zeros((dimension, dimension))

    for n in range(0, dimension):
        mat[n][n] = numpy.random.randint(a, b + 1)

    Q = rvs(loc, scale, dimension)
    M = Q.dot(mat).dot(Q.T)

    return M


def generate_b_vector(a, b, dimension):
    b_vector = numpy.random.rand(dimension, 1)

    for n in range(0, dimension):
        b_vector[n] = (b - a) * b_vector[n] + a

    return b_vector

def main():
    if len(sys.argv) < 7:
        print('Usage:')
        print("{} output_data_path loc scale a b dimension".format(sys.argv[0]))
        sys.exit(1)

    output_data_path = sys.argv[1]
    loc = float(sys.argv[2])
    scale = float(sys.argv[3])
    a = float(sys.argv[4])
    b = float(sys.argv[5])
    dimension = int(sys.argv[6])

    matrix = generate_a_matrix(loc, scale, a, b, dimension)
    b_vector = generate_b_vector(a, b, dimension)

    data = numpy.append(matrix, b_vector, axis=1)

    with open(output_data_path, 'w') as file:
        for row in range(0, dimension):
            file.writelines(' '.join([format_output(x) for x in data[row]]) + '\n')

    solution = numpy.linalg.solve(matrix, b_vector)

    with open(output_data_path + '.expected_solution', 'w') as file:
        for n in solution:
            file.writelines(format_output(n) + '\n')


if __name__ == '__main__':
    main()
