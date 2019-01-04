cmake_minimum_required(VERSION 3.13)
project(chebyshev_tests C)

set(CMAKE_C_STANDARD 11)
set(NUMPROC 4)

file(GLOB SOURCES
        ../src/*.h
        ../src/*.c
        )

find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

add_executable(basic_test basic.c ${SOURCES})
add_test(basic_test ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${NUMPROC} basic_test)
target_link_libraries(basic_test PRIVATE m ${MPI_C_LIBRARIES})

add_executable(various_s_parameter_test various_s_parameter.c ${SOURCES})
add_test(various_s_parameter_test ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${NUMPROC} various_s_parameter_test)
target_link_libraries(various_s_parameter_test PRIVATE m ${MPI_C_LIBRARIES})

file(GLOB TEST_FILES
        "sources/*"
        )

file(COPY ${TEST_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/sources)
