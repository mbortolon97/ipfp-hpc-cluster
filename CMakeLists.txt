cmake_minimum_required(VERSION 3.1)
project(DistributedIPFP)

find_package(MPI REQUIRED)

add_executable(distributedIPFP src/main-ipfp.c src/dense_distribution.c src/dense_matrix.c src/matrix_partition.c src/permutation.c src/sparse_matrix.c src/submatrix.c src/log.c)
# src/dense_distribution.h src/dense_matrix.h src/message_tag.h src/matrix_partition.h src/permutation.h src/sparse_matrix.h
target_include_directories(distributedIPFP PRIVATE ${MPI_CXX_INCLUDE_PATH})
target_compile_options(distributedIPFP PRIVATE ${MPI_CXX_COMPILE_FLAGS})
target_link_libraries(distributedIPFP ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS})