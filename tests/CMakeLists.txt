cmake_minimum_required(VERSION 3.13)
project(chebyshev_tests C)

set(CMAKE_C_STANDARD 11)

file(GLOB SOURCES
        ../src/*.h
        ../src/*.c
        )

add_executable(basic_tests basic_tests.c ${SOURCES})
add_test(basic_tests basic_tests)
