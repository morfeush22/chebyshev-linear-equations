cmake_minimum_required(VERSION 3.13)
project(chebyshev C)

set(CMAKE_C_STANDARD 11)

enable_testing()
add_subdirectory(tests)

add_executable(main main.c)
