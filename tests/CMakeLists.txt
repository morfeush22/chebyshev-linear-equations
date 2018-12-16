cmake_minimum_required(VERSION 3.13)
project(chebyshev_tests C)

set(CMAKE_C_STANDARD 11)

add_executable(basic_tests ../src/base.c basic_tests.c)
add_test(basic_tests basic_tests)
