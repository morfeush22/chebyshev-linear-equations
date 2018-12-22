cmake_minimum_required(VERSION 3.13)
project(chebyshev_tests C)

set(CMAKE_C_STANDARD 11)

file(GLOB SOURCES
        ../src/*.h
        ../src/*.c
        )

add_executable(basic_test basic.c ${SOURCES})
add_test(basic_test basic_test)

add_executable(various_s_parameter_test various_s_parameter.c ${SOURCES})
add_test(various_s_parameter_test various_s_parameter_test)
