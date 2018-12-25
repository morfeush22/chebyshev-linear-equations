cmake_minimum_required(VERSION 3.13)
project(chebyshev C)

set(CMAKE_C_STANDARD 11)

enable_testing()
add_subdirectory(tests)

file(GLOB SOURCES
        src/*.h
        src/*.c
        )

find_package(OpenMP)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
endif()

set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -DDEBUG_PRINTS")

add_executable(main main.c ${SOURCES})
target_link_libraries(main m)
