cmake_minimum_required(VERSION 3.1)
project(DistributedIPFP)

find_package(MPI REQUIRED)
add_executable(my_mpi_bin src1.cpp src2.cpp)
target_include_directories(my_mpi_bin PRIVATE ${MPI_CXX_INCLUDE_PATH} src1.h src2.h)
target_compile_options(my_mpi_bin PRIVATE ${MPI_CXX_COMPILE_FLAGS} my_compile_flags)
target_link_libraries(my_mpi_bin ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS} my_link_flags)