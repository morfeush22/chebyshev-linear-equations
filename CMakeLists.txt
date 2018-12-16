cmake_minimum_required(VERSION 3.13)
project(chebyshev C)

set(CMAKE_C_STANDARD 11)

enable_testing()
add_subdirectory(tests)

file(GLOB SOURCES
        src/*.h
        src/*.c
        )

set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -DDEBUG_PRINTS")

add_executable(main main.c ${SOURCES})
