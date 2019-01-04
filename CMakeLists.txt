cmake_minimum_required(VERSION 3.13)
project(chebyshev C)

set(CMAKE_C_STANDARD 11)

enable_testing()
add_subdirectory(tests)

file(GLOB SOURCES
        src/*.h
        src/*.c
        )

find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -DDEBUG_PRINTS")

add_executable(main main.c ${SOURCES})
target_link_libraries(main PRIVATE m ${MPI_C_LIBRARIES})
