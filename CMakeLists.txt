cmake_minimum_required(VERSION 3.1)
project(DistributedIPFP)

find_package(MPI REQUIRED)

add_executable(distributedIPFP src/main-ipfp.c src/dense_matrix.c src/matrix_partition.c src/permutation.c src/smallmatrix.c src/sparse_matrix.c)
target_include_directories(distributedIPFP PRIVATE ${MPI_CXX_INCLUDE_PATH} src/dense_matrix.h src/matrix_partition.h src/permutation.h src/smallmatrix.h src/sparse_matrix.h)
target_compile_options(distributedIPFP PRIVATE ${MPI_CXX_COMPILE_FLAGS})
target_link_libraries(distributedIPFP ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS})